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
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed,TimeoutError
import pandas as pd
from tqdm import tqdm
import traceback
import json
import os

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
    mode = 'ligand' # all, structure, ligand
) -> None:
    """Evaluate the structure."""
    if mode == 'ligand':
        result = subprocess.run(
            OST_COMPARE_LIGAND_STRUCTURE.format(
                model_file=str(pred),
                reference_file=str(reference),
                output_path=os.path.join(outdir),
                ),
                shell=True,  # noqa: S602
                check=False,
                executable=executable,
                capture_output=True,
            )
    elif mode == 'structure':
        result = subprocess.run(
            OST_COMPARE_STRUCTURE.format(
                model_file=str(pred),
                reference_file=str(reference),
                output_path=os.path.join(outdir),
                ),
                shell=True,  # noqa: S602
                check=False,
                executable=executable,
                capture_output=True,
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

    if not os.path.exists(prediction_path):
        print(f"prediction_path is None for {pdb_id} with seed {seed} and sample {sample}")
        return "prediction_path is None"
    
    native_path = os.path.join(ground_truth_path, f'{pdb_id}.cif')
    
    try:
        evaluate_structure(
            pred = prediction_path,
            reference = native_path,
            outdir = output_path,
            mode = mode
        )
    except BaseException as e:
        print(f"Error when calculating dockq for {pdb_id} with seed {seed} and sample {sample}")
        print(traceback.format_exc())
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

    if not os.path.exists(prediction_path):
        print(f"prediction_path is None for {pdb_id} with seed {seed} and sample {sample}")
        return "prediction_path is None"


    result = {
            **row,
        }

    try:
        if mode == 'ligand':
            rmsd,lddt_lp, lddt_pli = get_ligand_value(output_path,interface_chain_id_1,interface_chain_id_2)
            result.update({
                'rmsd': rmsd,
                'lddt-lp': lddt_lp,
                'lddt-pli': lddt_pli,
            })
        elif mode == 'structure':
            dockq_ost,irmsd,lrmsd,len_dockq,lddt,tm_score,gdt_ts,rmsd = get_structure_value(output_path,interface_chain_id_1,interface_chain_id_2)
            result.update({
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
