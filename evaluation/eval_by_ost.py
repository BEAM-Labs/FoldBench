import itertools
from functools import partial
import pandas as pd
from tqdm import tqdm
import json
import os
from multiprocessing import Pool
import traceback
import subprocess
import numpy as np
import signal
import uuid
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed,TimeoutError
import pandas as pd
from tqdm import tqdm
import traceback
import json
import os

OST_TIMEOUT_SEC = int(os.environ.get("FOLDBENCH_OST_TIMEOUT_SEC", "600"))

OST_COMPARE_LIGAND_STRUCTURE = r"""
ost compare-ligand-structures \
-m {model_file} \
-r {reference_file} \
--fault-tolerant \
--lddt-pli --rmsd \
-o {output_path}
"""

OST_COMPARE_STRUCTURE = r"""
    ost compare-structures \
        -m {model_file} \
        -r {reference_file} \
        -o {output_path} \
        --fault-tolerant \
        --min-pep-length 4 \
        --min-nuc-length 4 \
        --lddt --rigid-scores --tm-score --dockq \
"""

def _run_ost_command(command, executable, timeout_sec):
    """Run an OST command and terminate the whole process group on timeout."""
    process = subprocess.Popen(
        command,
        shell=True,  # noqa: S602
        executable=executable,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout_sec)
    except subprocess.TimeoutExpired as exc:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            stdout, stderr = process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = process.communicate()
        exc.stdout = stdout
        exc.stderr = stderr
        raise exc

    return subprocess.CompletedProcess(
        args=command,
        returncode=process.returncode,
        stdout=stdout,
        stderr=stderr,
    )

def _safe_json_value(value):
    if pd.isna(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value

def _is_valid_path(value):
    if not isinstance(value, (str, bytes, os.PathLike)):
        return False
    return bool(os.fspath(value))

def write_failed_detail(output_path, row, mode, reason, message=None, stdout=None, stderr=None):
    failed_detail = {
        "status": "failed",
        "failure_reason": reason,
        "failure_message": message,
        "mode": mode,
        "pdb_id": _safe_json_value(row.get("pdb_id")),
        "seed": _safe_json_value(row.get("seed")),
        "sample": _safe_json_value(row.get("sample")),
        "prediction_path": _safe_json_value(row.get("prediction_path")),
    }
    if stdout:
        failed_detail["stdout_tail"] = stdout[-4000:]
    if stderr:
        failed_detail["stderr_tail"] = stderr[-4000:]

    tmp_output_path = f"{output_path}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
    with open(tmp_output_path, "w") as f:
        json.dump(failed_detail, f, indent=2)
    os.replace(tmp_output_path, output_path)

def read_failed_detail(output_path):
    try:
        with open(output_path, "r") as f:
            data = json.load(f)
    except FileNotFoundError:
        return None
    except json.JSONDecodeError:
        return {
            "status": "failed",
            "failure_reason": "invalid_json",
        }
    if data.get("status") == "failed":
        return data
    return None

def get_structure_value(output_path,native_chain_id_1,native_chain_id_2):

    # Load the CSV file
    dockq = None
    irmsd = None
    lrmsd = None
    len_dockq = None
    lddt = None
    tm_score = None
    gdt_ts = None
    rmsd = None
    try:
        with open(output_path, 'r') as f:
            data = json.load(f)
        

        for i,interface in enumerate(data['dockq_interfaces']):
            if native_chain_id_1 in interface and native_chain_id_2 in interface:
                dockq = data['dockq'][i]
                irmsd = data['irmsd'][i]
                lrmsd = data['lrmsd'][i]
                # breakpoint()
        len_dockq = len(data['dockq'])
        lddt = data['lddt']
        tm_score = data['tm_score']
        gdt_ts = data['oligo_gdtts']
        rmsd = data['rmsd']
    except FileNotFoundError:
        print(f"FileNotFoundError: {output_path}")  # Return None if the file does not exist
    finally:
        return dockq,irmsd,lrmsd,len_dockq,lddt,tm_score,gdt_ts,rmsd

def get_ligand_value(output_path,native_chain_id_1,native_chain_id_2):

    # Load the CSV file
    rmsd = None
    lddt_lp = None
    lddt_pli = None
    try:
        with open(output_path, 'r') as f:
            data = json.load(f)
        
        for item in data['rmsd']['assigned_scores']:
            reference_ligand_name = item['reference_ligand']
            #  native_chain_id_2 is ligand by default
            if native_chain_id_2 == reference_ligand_name.split('.')[0]:
                rmsd = item['score']
                lddt_lp = item['lddt_lp']
        
        for item in data['lddt_pli']['assigned_scores']:
            reference_ligand_name = item['reference_ligand']
            #  native_chain_id_2 is ligand by default
            if native_chain_id_2  == reference_ligand_name.split('.')[0]:
                lddt_pli = item['score']
    
    except FileNotFoundError:
        print(f"FileNotFoundError: {output_path}")  # Return None if the file does not exist
    finally:
        return rmsd,lddt_lp,lddt_pli

def evaluate_structure(
    pred,
    reference,
    outdir: str,
    executable: str = "/bin/bash",
    mode = 'ligand', # all, structure, ligand
    timeout_sec = OST_TIMEOUT_SEC,
) -> None:
    """Evaluate the structure."""
    if mode == 'ligand':
        result = _run_ost_command(
            OST_COMPARE_LIGAND_STRUCTURE.format(
                model_file=str(pred),
                reference_file=str(reference),
                output_path=os.path.join(outdir),
                ),
                executable=executable,
                timeout_sec=timeout_sec,
            )
    elif mode == 'structure':
        result = _run_ost_command(
            OST_COMPARE_STRUCTURE.format(
                model_file=str(pred),
                reference_file=str(reference),
                output_path=os.path.join(outdir),
                ),
                executable=executable,
                timeout_sec=timeout_sec,
            )
    return result

def ost_evaluation(args):
    row, ground_truth_path, detail_path, mode= args

    pdb_id = row["pdb_id"]
    seed = row["seed"]
    sample = row["sample"]
    prediction_path = row["prediction_path"]

    output_path = f'{detail_path}/{pdb_id}_{seed}_{sample}_{mode}_ost.json'
    if os.path.exists(output_path):
        return "exist"

    if not _is_valid_path(prediction_path) or not os.path.exists(prediction_path):
        print(f"prediction_path is None for {pdb_id} with seed {seed} and sample {sample}")
        write_failed_detail(
            output_path,
            row,
            mode,
            reason="prediction_path_missing",
            message=f"prediction_path is missing or does not exist: {prediction_path}",
        )
        return "prediction_path is None"
    
    native_path = os.path.join(ground_truth_path, f'{pdb_id}.cif')
    
    try:
        result = evaluate_structure(
            pred = prediction_path,
            reference = native_path,
            outdir = output_path,
            mode = mode
        )
        if result.returncode != 0 and not os.path.exists(output_path):
            write_failed_detail(
                output_path,
                row,
                mode,
                reason="ost_nonzero_returncode",
                message=f"OST exited with return code {result.returncode}",
                stdout=result.stdout,
                stderr=result.stderr,
            )
            return "ost_nonzero_returncode"
    except subprocess.TimeoutExpired as e:
        print(f"OST timeout for {pdb_id} with seed {seed} and sample {sample} after {OST_TIMEOUT_SEC}s")
        write_failed_detail(
            output_path,
            row,
            mode,
            reason="timeout",
            message=f"OST command exceeded {OST_TIMEOUT_SEC}s",
            stdout=getattr(e, "stdout", None),
            stderr=getattr(e, "stderr", None),
        )
        return "timeout"
    except BaseException as e:
        print(f"Error when calculating dockq for {pdb_id} with seed {seed} and sample {sample}")
        print(traceback.format_exc())
        write_failed_detail(
            output_path,
            row,
            mode,
            reason="exception",
            message=traceback.format_exc(),
        )
        return "Error when calculating dockq"

def ost_get_result(args):
    
    row, ground_truth_path, detail_path, mode= args

    pdb_id = row["pdb_id"]
    if 'interface_chain_id_1' in row.keys() and 'interface_chain_id_2' in row.keys():
        interface_chain_id_1 = row["interface_chain_id_1"]
        interface_chain_id_2 = row["interface_chain_id_2"]
    elif 'native_chain_id_1' in row.keys() and 'native_chain_id_2' in row.keys():
        interface_chain_id_1 = row["native_chain_id_1"]
        interface_chain_id_2 = row["native_chain_id_2"]
    elif 'chain_id' in row.keys():
        interface_chain_id_1 = row["chain_id"]
        interface_chain_id_2 = row["chain_id"]
    seed = row["seed"]
    sample = row["sample"]
    prediction_path = row["prediction_path"]
    
    output_path = f'{detail_path}/{pdb_id}_{seed}_{sample}_{mode}_ost.json'

    if not _is_valid_path(prediction_path) or not os.path.exists(prediction_path):
        print(f"prediction_path is None for {pdb_id} with seed {seed} and sample {sample}")
        result = {
                **row,
            }
        result.update({
            'evaluation_status': 'failed',
            'evaluation_failure_reason': 'prediction_path_missing',
            'evaluation_failure_message': f"prediction_path is missing or does not exist: {prediction_path}",
        })
        if mode == 'ligand':
            result.update({
                'rmsd': None,
                'lddt-lp': None,
                'lddt-pli': None,
            })
        elif mode == 'structure':
            result.update({
                'dockq_score': None,
                'irmsd': None,
                'lrmsd': None,
                'len_dockq': None,
                'lddt': None,
                'tm_score': None,
                'gdt_ts': None,
                'rmsd': None,
            })
        return result


    result = {
            **row,
        }

    try:
        failed_detail = read_failed_detail(output_path)
        if failed_detail is not None:
            result.update({
                'evaluation_status': 'failed',
                'evaluation_failure_reason': failed_detail.get('failure_reason'),
                'evaluation_failure_message': failed_detail.get('failure_message'),
            })
            if mode == 'ligand':
                result.update({
                    'rmsd': None,
                    'lddt-lp': None,
                    'lddt-pli': None,
                })
            elif mode == 'structure':
                result.update({
                    'dockq_score': None,
                    'irmsd': None,
                    'lrmsd': None,
                    'len_dockq': None,
                    'lddt': None,
                    'tm_score': None,
                    'gdt_ts': None,
                    'rmsd': None,
                })
            return result

        if mode == 'ligand':
            rmsd,lddt_lp, lddt_pli = get_ligand_value(output_path,interface_chain_id_1,interface_chain_id_2)
            result.update({
                'evaluation_status': 'success',
                'evaluation_failure_reason': None,
                'evaluation_failure_message': None,
                'rmsd': rmsd,
                'lddt-lp': lddt_lp,
                'lddt-pli': lddt_pli,
            })
        elif mode == 'structure':
            dockq_ost,irmsd,lrmsd,len_dockq,lddt,tm_score,gdt_ts,rmsd = get_structure_value(output_path,interface_chain_id_1,interface_chain_id_2)
            result.update({
                'evaluation_status': 'success',
                'evaluation_failure_reason': None,
                'evaluation_failure_message': None,
                'dockq_score': dockq_ost,
                'irmsd':irmsd,
                'lrmsd':lrmsd,
                'len_dockq':len_dockq,
                'lddt':lddt,
                'tm_score':tm_score,
                'gdt_ts':gdt_ts,
                'rmsd':rmsd
            })

        return result
    except BaseException as e:
        print(f"Error when calculating dockq for {pdb_id} with seed {seed} and sample {sample}")
        print(traceback.format_exc())
        return None



def eval_by_ost(target_df,target_type,evaluation_dir,ground_truth_dir,max_workers = 64):

    detail_path = os.path.join(evaluation_dir, 'detail')
    if not os.path.exists(detail_path):
        os.makedirs(detail_path)

    mode = ''
    
    if target_type in ["interface_protein_protein","interface_antibody_antigen","interface_protein_peptide","interface_protein_dna","interface_protein_rna","monomer_dna","monomer_rna","monomer_protein"]:
        mode = "structure"
    elif target_type == "interface_protein_ligand":
        mode = "ligand"
    
    tasks = []
    for index, row in target_df.iterrows():
        tasks.append((
                row,
                ground_truth_dir,
                detail_path,
                mode
            ))
    
    results_ost = []
    # evaluation by ost
    with ThreadPoolExecutor(max_workers=max_workers) as executor:

            future_to_task = {executor.submit(ost_evaluation , task): task for task in tasks}
            for future in tqdm(as_completed(future_to_task), total=len(tasks)):
                try:
                    result = future.result(timeout=20)
                    if result is not None:
                        results_ost.append(result)
                except TimeoutError:
                    print("this took too long...")
                    task = future_to_task[future]
                    future.cancel()
                except Exception as e:
                    task = future_to_task[future]
                    print(f"Error occurred for task: {task}")
                    print(traceback.format_exc())
                    future.cancel()
    
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {executor.submit(ost_get_result, task): task for task in tasks}

        for future in tqdm(as_completed(future_to_task), total=len(tasks)):
            try:
                result = future.result(timeout=20)
                if result is not None:  
                    results.append(result)
            except TimeoutError:
                print("this took too long...")
                task = future_to_task[future]
                future.cancel()#
            except Exception as e:
                task = future_to_task[future]
                print(f"Error occurred for task: {task}")
                print(traceback.format_exc())
                future.cancel()

    print(f"Total results for {target_type}: {len(results)}")
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(evaluation_dir,'raw',f"{target_type}_ost.csv"), index=False)
