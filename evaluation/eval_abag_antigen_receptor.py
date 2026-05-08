"""Antigen-centered antibody-antigen scoring.

The default ABAG evaluation continues to use ``ost compare-structures`` via
``eval_by_ost``. This module provides an explicit alternative scorer for the
case where the antigen should be treated as the fixed receptor and the antibody
chain should be scored as the moving ligand.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any
import json
import math
import os
import traceback
import uuid

import numpy as np
import pandas as pd
from tqdm import tqdm


FNAT_THRESHOLD = 5.0
INTERFACE_THRESHOLD = 10.0
BACKBONE_ATOMS_SEL = "aname=N,CA,C,O"
RAW_SUFFIX = "antigen_receptor_ost"


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj: Any) -> Any:
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def _safe_json_value(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value


def _safe_name(value: Any) -> str:
    text = str(value)
    for old, new in [("/", "_"), (" ", "_"), (",", "-"), (";", "-"), ("|", "-")]:
        text = text.replace(old, new)
    return text


def _is_valid_path(value: Any) -> bool:
    if not isinstance(value, (str, bytes, os.PathLike)):
        return False
    return bool(os.fspath(value))


def _quote_cnames(chains_str: str) -> str:
    parts = [part.strip() for part in str(chains_str).split(",") if part.strip()]
    return ",".join(f"'{part}'" for part in parts)


def _split_chain_names(chains_str: str) -> list[str]:
    return [part.strip() for part in str(chains_str).split(",") if part.strip()]


def _chain_residue_names(chain: Any) -> list[str]:
    return [str(residue.name) for residue in chain.residues]


def _chain_similarity(ref_chain: Any, model_chain: Any) -> float:
    ref_names = _chain_residue_names(ref_chain)
    model_names = _chain_residue_names(model_chain)
    if not ref_names or not model_names:
        return 0.0
    return float(SequenceMatcher(None, ref_names, model_names, autojunk=False).ratio())


def _infer_model_chain_map(ref: Any, model: Any) -> dict[str, str]:
    """Map reference chain names onto model chain names.

    AF3 may rename repeated assembly chains, e.g. reference ``A,B,A-2,B-2`` can
    become model ``A,B,F,G``. Targets use reference names, so this mapping is
    required before selecting model chains.
    """

    ref_chains = list(ref.chains)
    model_chains = list(model.chains)
    model_by_name = {chain.name: chain for chain in model_chains}
    used_model_names: set[str] = set()
    mapping: dict[str, str] = {}

    for ref_chain in ref_chains:
        model_chain = model_by_name.get(ref_chain.name)
        if model_chain is None:
            continue
        if _chain_similarity(ref_chain, model_chain) >= 0.60:
            mapping[ref_chain.name] = model_chain.name
            used_model_names.add(model_chain.name)

    for ref_chain in ref_chains:
        if ref_chain.name in mapping:
            continue
        best_name = None
        best_score = -1.0
        for model_chain in model_chains:
            if model_chain.name in used_model_names:
                continue
            score = _chain_similarity(ref_chain, model_chain)
            if score > best_score:
                best_score = score
                best_name = model_chain.name
        if best_name is not None and best_score >= 0.35:
            mapping[ref_chain.name] = best_name
            used_model_names.add(best_name)

    return mapping


def _translate_chains(chains_str: str, chain_map: dict[str, str]) -> str:
    translated = []
    missing = []
    for chain_name in _split_chain_names(chains_str):
        mapped = chain_map.get(chain_name)
        if mapped is None:
            missing.append(chain_name)
        else:
            translated.append(mapped)
    if missing:
        raise ValueError(f"unable to map reference chain(s) to model chain(s): {','.join(missing)}")
    return ",".join(translated)


def _match_views(model_view: Any, ref_view: Any) -> tuple[Any, Any]:
    from ost.mol import alg

    matched_res = alg.MatchResidueByNum(model_view, ref_view)
    if isinstance(matched_res, tuple):
        return matched_res[0], matched_res[1]
    return matched_res.view1, matched_res.view2


def _min_residue_distance_sq(res1: Any, res2: Any) -> float:
    from ost import geom

    min_dist_sq = float("inf")
    for atom_1 in res1.atoms:
        for atom_2 in res2.atoms:
            dist_sq = geom.Length2(atom_1.pos - atom_2.pos)
            if dist_sq < min_dist_sq:
                min_dist_sq = dist_sq
    return min_dist_sq


def _svd_rmsd(coords1: list[list[float]], coords2: list[list[float]]) -> float:
    p = np.array(coords1, dtype=np.float64)
    q = np.array(coords2, dtype=np.float64)
    p_center = p.mean(axis=0)
    q_center = q.mean(axis=0)
    p_c = p - p_center
    q_c = q - q_center
    h_mat = p_c.T @ q_c
    u_mat, _s_mat, vt_mat = np.linalg.svd(h_mat)
    det_val = np.linalg.det(vt_mat.T @ u_mat.T)
    sign_mat = np.diag([1.0, 1.0, np.sign(det_val)])
    rot = vt_mat.T @ sign_mat @ u_mat.T
    rotated = p_c @ rot.T
    diff = rotated - q_c
    return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))


def _dockq_formula(fnat: float, irmsd: float, lrmsd: float) -> float:
    return (
        float(fnat)
        + 1.0 / (1.0 + (irmsd / 1.5) ** 2)
        + 1.0 / (1.0 + (lrmsd / 8.5) ** 2)
    ) / 3.0


def manual_antigen_centered_scores(
    model_path: str,
    ref_path: str,
    antigen_chains: str,
    antibody_chains: str,
) -> dict[str, Any]:
    from ost import geom, io
    from ost.mol import alg

    model = io.LoadMMCIF(model_path).Select("peptide=True")
    ref = io.LoadMMCIF(ref_path).Select("peptide=True")

    chain_map = _infer_model_chain_map(ref, model)
    model_antigen_chains = _translate_chains(antigen_chains, chain_map)
    model_antibody_chains = _translate_chains(antibody_chains, chain_map)

    ref_rec_sel = f"cname={_quote_cnames(antigen_chains)}"
    ref_lig_sel = f"cname={_quote_cnames(antibody_chains)}"
    model_rec_sel = f"cname={_quote_cnames(model_antigen_chains)}"
    model_lig_sel = f"cname={_quote_cnames(model_antibody_chains)}"

    m_rec_match, r_rec_match = _match_views(model.Select(model_rec_sel), ref.Select(ref_rec_sel))
    m_lig_match, r_lig_match = _match_views(model.Select(model_lig_sel), ref.Select(ref_lig_sel))

    r_rec_residues = list(r_rec_match.residues)
    r_lig_residues = list(r_lig_match.residues)
    m_rec_residues = list(m_rec_match.residues)
    m_lig_residues = list(m_lig_match.residues)

    if not r_rec_residues or not r_lig_residues:
        raise ValueError(
            "empty receptor or ligand selection: "
            f"antigen={antigen_chains}, antibody={antibody_chains}"
        )
    if len(r_rec_residues) != len(m_rec_residues) or len(r_lig_residues) != len(m_lig_residues):
        raise ValueError(
            "matched residue count mismatch: "
            f"antigen ref/model={len(r_rec_residues)}/{len(m_rec_residues)}, "
            f"antibody ref/model={len(r_lig_residues)}/{len(m_lig_residues)}"
        )

    native_contacts: set[tuple[int, int]] = set()
    rec_interface_idx: set[int] = set()
    lig_interface_idx: set[int] = set()
    fnat_thresh_sq = FNAT_THRESHOLD**2
    iface_thresh_sq = INTERFACE_THRESHOLD**2

    for i_res, r_rec_res in enumerate(r_rec_residues):
        for j_res, r_lig_res in enumerate(r_lig_residues):
            dist_sq = _min_residue_distance_sq(r_rec_res, r_lig_res)
            if dist_sq < fnat_thresh_sq:
                native_contacts.add((i_res, j_res))
            if dist_sq < iface_thresh_sq:
                rec_interface_idx.add(i_res)
                lig_interface_idx.add(j_res)

    nat_total = len(native_contacts)
    if nat_total == 0:
        return {
            "dockq_score": 0.0,
            "fnat": 0.0,
            "irmsd": None,
            "lrmsd": None,
            "nat_correct": 0,
            "nat_total": 0,
            "quality": "Incorrect",
            "model_antigen_chain": model_antigen_chains,
            "model_antibody_chain": model_antibody_chains,
            "chain_map": chain_map,
        }

    nat_correct = 0
    for i_res, j_res in native_contacts:
        dist_sq = _min_residue_distance_sq(m_rec_residues[i_res], m_lig_residues[j_res])
        if dist_sq < fnat_thresh_sq:
            nat_correct += 1
    fnat = nat_correct / nat_total

    bb_names = ["N", "CA", "C", "O"]
    model_iface_coords: list[list[float]] = []
    ref_iface_coords: list[list[float]] = []
    for i_res in sorted(rec_interface_idx):
        for atom_name in bb_names:
            model_atom = m_rec_residues[i_res].FindAtom(atom_name)
            ref_atom = r_rec_residues[i_res].FindAtom(atom_name)
            if model_atom.IsValid() and ref_atom.IsValid():
                model_iface_coords.append([model_atom.pos.x, model_atom.pos.y, model_atom.pos.z])
                ref_iface_coords.append([ref_atom.pos.x, ref_atom.pos.y, ref_atom.pos.z])
    for j_res in sorted(lig_interface_idx):
        for atom_name in bb_names:
            model_atom = m_lig_residues[j_res].FindAtom(atom_name)
            ref_atom = r_lig_residues[j_res].FindAtom(atom_name)
            if model_atom.IsValid() and ref_atom.IsValid():
                model_iface_coords.append([model_atom.pos.x, model_atom.pos.y, model_atom.pos.z])
                ref_iface_coords.append([ref_atom.pos.x, ref_atom.pos.y, ref_atom.pos.z])

    irmsd = _svd_rmsd(model_iface_coords, ref_iface_coords) if model_iface_coords else float("inf")

    m_rec_bb = m_rec_match.Select(BACKBONE_ATOMS_SEL)
    r_rec_bb = r_rec_match.Select(BACKBONE_ATOMS_SEL)
    alg.SuperposeSVD(m_rec_bb, r_rec_bb)

    m_lig_bb = m_lig_match.Select(BACKBONE_ATOMS_SEL)
    r_lig_bb = r_lig_match.Select(BACKBONE_ATOMS_SEL)
    sum_sq_dist = 0.0
    atom_count = 0
    for atom_model, atom_ref in zip(m_lig_bb.atoms, r_lig_bb.atoms):
        sum_sq_dist += geom.Length2(atom_model.pos - atom_ref.pos)
        atom_count += 1
    lrmsd = math.sqrt(sum_sq_dist / atom_count) if atom_count else float("inf")

    dockq = _dockq_formula(fnat, irmsd, lrmsd)
    if dockq >= 0.80:
        quality = "High"
    elif dockq >= 0.49:
        quality = "Medium"
    elif dockq >= 0.23:
        quality = "Acceptable"
    else:
        quality = "Incorrect"

    return {
        "dockq_score": round(float(dockq), 3),
        "fnat": round(float(fnat), 3),
        "nat_correct": int(nat_correct),
        "nat_total": int(nat_total),
        "irmsd": round(float(irmsd), 3) if math.isfinite(irmsd) else None,
        "lrmsd": round(float(lrmsd), 3) if math.isfinite(lrmsd) else None,
        "quality": quality,
        "model_antigen_chain": model_antigen_chains,
        "model_antibody_chain": model_antibody_chains,
        "chain_map": chain_map,
    }


def _write_detail(output_path: Path, payload: dict[str, Any]) -> None:
    tmp_output_path = output_path.with_suffix(f"{output_path.suffix}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    with open(tmp_output_path, "w") as handle:
        json.dump(payload, handle, indent=2, cls=NumpyEncoder)
    os.replace(tmp_output_path, output_path)


def _failed_result(row: dict[str, Any], reason: str, message: str | None) -> dict[str, Any]:
    result = {**row}
    result.update(
        {
            "evaluation_status": "failed",
            "evaluation_failure_reason": reason,
            "evaluation_failure_message": message,
            "dockq_score": None,
            "fnat": None,
            "nat_correct": None,
            "nat_total": None,
            "irmsd": None,
            "lrmsd": None,
            "quality": None,
            "model_antigen_chain": None,
            "model_antibody_chain": None,
            "chain_map": None,
        }
    )
    return result


def _evaluate_one(row: dict[str, Any], ground_truth_dir: str, detail_path: str) -> dict[str, Any]:
    pdb_id = row["pdb_id"]
    seed = row.get("seed")
    sample = row.get("sample")
    prediction_path = row.get("prediction_path")
    antibody_chain = str(row.get("antibody_chain", ""))
    antigen_chain = str(row.get("antigen_chain", ""))

    detail_name = (
        f"{_safe_name(pdb_id)}_ab{_safe_name(antibody_chain)}_ag{_safe_name(antigen_chain)}"
        f"_seed{_safe_name(seed)}_sample{_safe_name(sample)}_{RAW_SUFFIX}.json"
    )
    detail_file = Path(detail_path) / detail_name

    if not antibody_chain or antibody_chain == "nan" or not antigen_chain or antigen_chain == "nan":
        result = _failed_result(row, "missing_abag_chain_columns", "antibody_chain/antigen_chain is required")
        _write_detail(detail_file, result)
        return result

    native_path = Path(ground_truth_dir) / f"{pdb_id}.cif"
    if not native_path.exists():
        result = _failed_result(row, "native_file_missing", str(native_path))
        _write_detail(detail_file, result)
        return result

    if not _is_valid_path(prediction_path) or not Path(prediction_path).exists():
        result = _failed_result(
            row,
            "prediction_path_missing",
            f"prediction_path is missing or does not exist: {prediction_path}",
        )
        _write_detail(detail_file, result)
        return result

    try:
        scores = manual_antigen_centered_scores(
            model_path=str(prediction_path),
            ref_path=str(native_path),
            antigen_chains=antigen_chain,
            antibody_chains=antibody_chain,
        )
        result = {**row}
        result.update(
            {
                "evaluation_status": "success",
                "evaluation_failure_reason": None,
                "evaluation_failure_message": None,
                "dockq_score": scores["dockq_score"],
                "fnat": scores["fnat"],
                "nat_correct": scores["nat_correct"],
                "nat_total": scores["nat_total"],
                "irmsd": scores["irmsd"],
                "lrmsd": scores["lrmsd"],
                "quality": scores["quality"],
                "model_antigen_chain": scores["model_antigen_chain"],
                "model_antibody_chain": scores["model_antibody_chain"],
                "chain_map": json.dumps(scores["chain_map"], sort_keys=True),
            }
        )
        _write_detail(detail_file, {**result, "scores": scores})
        return result
    except Exception as exc:
        result = _failed_result(row, exc.__class__.__name__, str(exc))
        result["traceback"] = traceback.format_exc()
        _write_detail(detail_file, result)
        return result


def eval_by_abag_antigen_receptor(target_df, target_type, evaluation_dir, ground_truth_dir, max_workers=64):
    if target_type != "interface_antibody_antigen":
        raise ValueError("Antigen-receptor scoring is only valid for interface_antibody_antigen")
    if "antibody_chain" not in target_df.columns or "antigen_chain" not in target_df.columns:
        raise ValueError("Antigen-receptor ABAG scoring requires antibody_chain and antigen_chain columns")

    detail_path = os.path.join(evaluation_dir, "detail")
    raw_path = os.path.join(evaluation_dir, "raw")
    os.makedirs(detail_path, exist_ok=True)
    os.makedirs(raw_path, exist_ok=True)

    rows = [{key: _safe_json_value(value) for key, value in row.to_dict().items()} for _, row in target_df.iterrows()]
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(_evaluate_one, row, ground_truth_dir, detail_path)
            for row in rows
        ]
        for future in tqdm(as_completed(futures), total=len(futures)):
            results.append(future.result())

    print(f"Total antigen-receptor ABAG results for {target_type}: {len(results)}")
    df = pd.DataFrame(results)
    output_path = os.path.join(raw_path, f"{target_type}_{RAW_SUFFIX}.csv")
    df.to_csv(output_path, index=False)
