"""
postprocess for algorithm. 
1. It is required to generate prediction_summary.csv in ./outputs/evaluation/{algorithm} for our benchmark. The keys should include pdb_id, seed, sample, ranking_score, and prediction_path.
2. Convert predictions from the algorithm output format to a format supported by OpenStructure and DockQv2. If the sample test succeeds, you can obtain the scores for each target in ./outputs/evaluation/{algorithm}.
You can use PostProcess.postprocess() to perform postprocessing.  
"""

import argparse
import os
import re
import numpy as np
from tqdm import tqdm
import json
import biotite.structure.io.pdbx as pdbx
import glob
import multiprocessing as mp    
import pandas as pd


class PostProcess():
    def __init__(self):
        pass
    
    def process_file(self,cif_paths):
        pdb_id,file_path,new_file_path,seed,sample = cif_paths

        if not os.path.exists(file_path):
            return None

        dict_entity_id = {}
        
        cif_file = pdbx.CIFFile.read(file_path)
        block = cif_file.block
        atom_site = block.get("atom_site")

        atom_site["occupancy"] = pdbx.CIFColumn(pdbx.CIFData(["1" for _ in range(len(atom_site['group_PDB']))]))
        atom_site['B_iso_or_equiv'] = pdbx.CIFColumn(pdbx.CIFData(["0" for _ in range(len(atom_site['group_PDB']))]))

        label_entity_ids = atom_site['label_entity_id'].as_array().tolist()
        group_pdbs = atom_site['group_PDB'].as_array().tolist()
        for entity_id, group_pdb in zip(label_entity_ids, group_pdbs):
            if entity_id not in dict_entity_id:
                if group_pdb == 'ATOM':
                    dict_entity_id[entity_id] = ('polymer', group_pdb)
                else:
                    dict_entity_id[entity_id] = ('non-polymer', group_pdb)
            else:  # modification
                if dict_entity_id[entity_id][1] != group_pdb:
                    dict_entity_id[entity_id] = ('polymer', group_pdb)
        
        block['entity'] = pdbx.CIFCategory({"id": list(dict_entity_id.keys()), 
                                            "type": [dict_entity_id[key][0] for key in dict_entity_id.keys()]})
        
        cif_file.write(new_file_path)
        
        return pdb_id,new_file_path,seed,sample

    def postprocess(self,input_dir,prediction_dir,evaluation_dir):
        
        with open(os.path.join(input_dir, "inputs.json"), "r") as f:
            input_data = json.load(f)

        seeds = ['42','66','101','2024','8888']
        samples = [0,1,2,3,4]
        cif_paths = []

        for input_data in input_data:
            pdb_id = input_data['name']
            for seed in seeds:
                for sample in samples:
                    cif_path = f"{prediction_dir}/{pdb_id}/seed_{seed}/predictions/{pdb_id}_seed_{seed}_sample_{sample}.cif"
                    # print(f"Processing {cif_path}")
                    
                    if os.path.exists(cif_path):
                        cif_new_path =  os.path.join(os.path.dirname(cif_path),
                                f'{os.path.splitext(os.path.basename(cif_path))[0]}_postprocessed.cif')
                        cif_paths.append((pdb_id,cif_path,cif_new_path,seed,sample))
        
        print(f"Processing {len(cif_paths)} files")
        num_cores = mp.cpu_count()
        # Use approximately 75% of CPU cores to avoid system overload
        num_processes = max(1, int(num_cores * 0.8))
        print(f"Will use {num_processes} processes for parallel processing")
        
        # Create process pool
        with mp.Pool(processes=num_processes) as pool:
            results = list(tqdm(
                pool.imap(self.process_file, cif_paths), 
                total=len(cif_paths),
                desc="Processing progress"
            ))
        
        data = []
        for result in results:
            if result is None:
                continue
            tmp = {}
            pdb_id,new_file_path,seed,sample = result
            
            # get ranking score
            confidence_path = f"{prediction_dir}/{pdb_id}/seed_{seed}/predictions/{pdb_id}_seed_{seed}_summary_confidence_sample_{sample}.json"
            with open(confidence_path, "r") as f:
                confidence = json.load(f)
            
            tmp['pdb_id'] = pdb_id
            tmp['seed'] = seed
            tmp['sample'] = sample
            tmp['ranking_score'] = confidence.get('ranking_score', 0)
            tmp['prediction_path'] = new_file_path
            
           
            data.append(tmp)

        df = pd.DataFrame(data)
        df.to_csv(os.path.join(evaluation_dir, f'prediction_summary.csv'), index=False)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--input_dir", required=True, help="The path to the algorithm input file."
)
parser.add_argument(
    "--prediction_dir", required=True, help="The path to the algorithm predictions file."
)
parser.add_argument(
    "--evaluation_dir", required=True, help="The dir with the evaluation files.",
)
args = parser.parse_args()


postprocess = PostProcess()

postprocess.postprocess(args.input_dir, args.prediction_dir,args.evaluation_dir)