

from evaluation import eval_by_dockqv2,eval_by_ost
import pandas as pd
import argparse
import os



parser = argparse.ArgumentParser()
parser.add_argument(
    "--targets_dir", required=False, default='./targets', help="The dir with the targets files."
)
parser.add_argument(
    "--prediction_dir", required=True, help="The path to the algorithm predictions file."
)
parser.add_argument(
    "--evaluation_dir", required=True,default='./outputs/evaluation', help="The dir with the evaluation files.",
)
parser.add_argument(
    "--ground_truth_dir", required=False, default='./ground_truths', help="The dir with the ground truth files.",
)
args = parser.parse_args()

# target_types = ["interface_protein_ligand","interface_protein_protein","interface_protein_peptide","interface_antibody_antigen","interface_protein_dna", "interface_protein_rna","monomer_dna","monomer_rna","monomer_protein"]
target_types = ["interface_protein_ligand","interface_antibody_antigen","interface_protein_dna","monomer_protein"]

prediction_summary_path = f'{args.evaluation_dir}/prediction_summary.csv'
prediction_summary_df = pd.read_csv(prediction_summary_path)

# caculation
for target_type in target_types:
    target_df_path = f'{args.targets_dir}/{target_type}.csv'
    if not os.path.exists(target_df_path):
        print(f"target_df_path is not exists for {target_type}")
        continue
    target_df = pd.read_csv(target_df_path)

    target_df = pd.merge(prediction_summary_df,target_df, on='pdb_id', how='left')

    if target_type in  ["interface_protein_protein","interface_antibody_antigen","interface_protein_peptide","interface_protein_ligand","interface_protein_dna","interface_protein_rna","monomer_dna","monomer_rna","monomer_protein"]:
        eval_by_ost(target_df,target_type,args.evaluation_dir,args.ground_truth_dir)
    
        if target_type in  ["interface_protein_dna","interface_protein_rna"]:
            eval_by_dockqv2(target_df,target_type,args.evaluation_dir,args.ground_truth_dir)