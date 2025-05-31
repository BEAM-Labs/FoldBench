import sys
from .DockQv2.DockQ import (
    load_PDB,
    run_on_all_native_interfaces,
    group_chains,
    get_all_chain_maps,
    count_chain_combinations,
    format_mapping
)
import itertools
from functools import partial
import pandas as pd
from tqdm import tqdm
import json
import os
from multiprocessing import Pool
import traceback
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Non-standard amino acids
    'MSE', 'SEC', 'PYL'
}

NUCLEOTIDES = {
    # DNA
    'DA', 'DT', 'DG', 'DC', 'DI', 'DU',
    # RNA
    'A', 'U', 'G', 'C', 'I'
}

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NumpyEncoder, self).default(obj)

def determine_chain_type(chain_id, residues):
    """
    Determine the type of a chain based on its residue composition.
    
    Parameters:
    - chain_id (str): Chain identifier
    - residues (list): List of residue names
    
    Returns:
    - str: Chain type ('protein', 'na', or 'unk')
    """
    if not residues:
        return 'unk'
        
    # Calculate the number of protein and nucleotide residues
    protein_count = sum(1 for res in residues if res in AMINO_ACIDS)
    na_count = sum(1 for res in residues if res in NUCLEOTIDES)
    
    total_residues = len(residues)
    
    # Set threshold (if more than 90% of residues are of a type, it is considered that type)
    threshold = 0.8
    
    # Determine the type
    if protein_count / total_residues >= threshold:
        return 'protein'
    elif na_count / total_residues >= threshold:
        return 'na'
    else:
        return 'unk'

# Add chain type to the structure
def reformat_type(structure):
    for chain_id, chain_value in structure.child_dict.items():
        residues=[]
        type='ligand' if chain_value.is_het else None
        if type is None:
            for res_id, res_value in chain_value.child_dict.items():
                residues.append(str(res_value.resname).upper())
            type=determine_chain_type(chain_id, residues)
        structure.child_dict[chain_id].type=type
        for idx, res in enumerate(structure.child_list):
                if res.id == chain_id:
                    structure.child_list[idx].type=type
    return structure


# Original DockQ treats modified residues as ligand, we need to reformat the structure to keep them as polymer
# NOTE: DockQv2.parsers.MMCIFParser has been modified to obtain polymer chain IDs from the structure
def reformat_het(structure):
    for chain_id, chain_value in structure.child_dict.items():
        is_polymer = False
        for res_id, res_value in chain_value.child_dict.items():
            if res_id[0] == ' ':
                is_polymer = True
                break
        if is_polymer:
            structure.child_dict[chain_id].is_het = None
            for idx, res in enumerate(structure.child_list):
                if res.id == chain_id:
                    structure.child_list[idx].is_het = None
                    break

    return structure



def dockq(model_path, native_path, model_chains=None, native_chains=None, small_molecule=False, allowed_mismatches=0):
    """
    Calculate the DockQ scores for a predicted structure.

    Parameters:
    - model_path (str): The path to the model (pred) PDB file.
    - native_path (str): The path to the native (ground truth) PDB file.
    - model_chains (list): A list of chain IDs in the model structure to consider. If None, all chains will be considered.
    - native_chains (list): A list of chain IDs in the native structure to consider. If None, all chains will be considered.
    - small_molecule (bool): Whether the structure contains a small molecule ligand. Default is False.
    - allowed_mismatches (int): The maximum number of allowed mismatches between model and native chains. Default is 0.
    """

    initial_mapping = {}

    model_structure = load_PDB(
        model_path, small_molecule=small_molecule
    )

    native_structure = load_PDB(
        native_path, small_molecule=small_molecule
    )

    native_structure = reformat_het(native_structure)
    model_structure = reformat_het(model_structure)
    model_structure = reformat_type(model_structure)
    native_structure = reformat_type(native_structure)

    model_chains = [
        c.id for c in model_structure] if model_chains is None else model_chains
    native_chains = [
        c.id for c in native_structure] if native_chains is None else native_chains

    # permute chains and run on a for loop
    best_dockq = -1
    best_result = None
    best_mapping = None

    model_chains_to_combo = [
        mc for mc in model_chains if mc not in initial_mapping.values()
    ]
    native_chains_to_combo = [
        nc for nc in native_chains if nc not in initial_mapping.keys()
    ]

    chain_clusters, reverse_map = group_chains(
        model_structure,
        native_structure,
        model_chains_to_combo,
        native_chains_to_combo,
        allowed_mismatches=allowed_mismatches
    )
    chain_maps = get_all_chain_maps(
        chain_clusters,
        initial_mapping,
        reverse_map,
        model_chains_to_combo,
        native_chains_to_combo,
    )

    num_chain_combinations = count_chain_combinations(chain_clusters)
    # copy iterator to use later
    chain_maps, chain_maps_ = itertools.tee(chain_maps)

    run_chain_map = partial(
        run_on_all_native_interfaces,
        model_structure,
        native_structure
    )

    if num_chain_combinations > 1:
        cpus = 1
        chunk_size = 1

        result_this_mappings = [run_chain_map(
            chain_map) for chain_map in chain_maps]

        for chain_map, (result_this_mapping, total_dockq) in zip(
            chain_maps_, result_this_mappings
        ):
            if total_dockq > best_dockq:
                best_dockq = total_dockq
                best_result = result_this_mapping
                best_mapping = chain_map

    else:
        best_mapping = next(chain_maps)
        best_result, best_dockq = run_chain_map(best_mapping)

    info = dict()
    info["model"] = model_path.split("/")[-1]
    info["native"] = native_path.split("/")[-1]
    info["best_dockq"] = best_dockq
    info["best_result"] = best_result
    info["GlobalDockQ"] = best_dockq / len(best_result)
    info["best_mapping"] = best_mapping

    return info


def process_single_case(args):
    
    row, ground_truth_path, detail_path, mode= args

    pdb_id = row["pdb_id"]
    interface_chain_id_1 = row["interface_chain_id_1"]
    interface_chain_id_2 = row["interface_chain_id_2"]
    seed = row["seed"]
    sample = row["sample"]
    prediction_path = row["prediction_path"]

    output_path = f'{detail_path}/{pdb_id}_{seed}_{sample}_{interface_chain_id_1}_{interface_chain_id_2}_{mode}_dockqv2.json'

    if not os.path.exists(prediction_path):
        print(f"prediction_path is None for {pdb_id} with seed {seed} and sample {sample}")
        return "prediction_path is None"
    
    native_path = os.path.join(ground_truth_path, f'{pdb_id}.cif')

    result = {
        **row,
    }
    
    if mode == 'ligand':
        small_molecule = True
    else:
        small_molecule = False

    try:
        info = dockq(
            model_path=prediction_path,
            native_path=native_path,
            native_chains=[interface_chain_id_1, interface_chain_id_2],
            small_molecule=small_molecule,
            allowed_mismatches=4
        )

        if info is None:
            return None
        else:
            json.dump(info, open(output_path, 'w'), cls=NumpyEncoder)

        
        key = list(info['best_result'].keys())[0]
        best_result = info["best_result"][key]
        
        result.update({
            'lrmsd': best_result["LRMSD"],
            'irmsd': best_result["iRMSD"],
            'dockq_score': best_result["DockQ"]
        })
        
        return result
        
    except BaseException as e:
        print(f"Error when calculating dockq for {pdb_id} with seed {seed} and sample {sample}")
        print(traceback.format_exc())
        return None







def eval_by_dockqv2(target_df,interface_type,evaluation_dir,ground_truth_dir,max_workers = 32):
    
    exported_path = evaluation_dir
    detail_path = os.path.join(exported_path, 'detail')
    if not os.path.exists(detail_path):
        os.makedirs(detail_path)
    mode = ''
    if interface_type in ["interface_protein_dna","interface_protein_rna"]:
        mode = "structure"
    elif interface_type == "ligand":
        mode = "ligand"
    
    tasks = []
    for index, row in target_df.iterrows():
        tasks.append((
                row,
                ground_truth_dir,
                detail_path,
                mode
            ))
    
    
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:

        future_to_task = {executor.submit(process_single_case, task): task for task in tasks}

        for future in tqdm(as_completed(future_to_task), total=len(tasks)):
            try:
                result = future.result(timeout=20)
                if result is not None:  
                    results.append(result)
            except TimeoutError:
                print("this took too long...")
                task = future_to_task[future]
                future.cancel()
            except Exception as e:
                task = future_to_task[future]
                print(f"Error occurred for task: {task}")
                print(traceback.format_exc())
                future.cancel()

    print(f"Total results for {interface_type}: {len(results)}")
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(evaluation_dir,'raw',f"{interface_type}_dockqv2.csv"), index=False)