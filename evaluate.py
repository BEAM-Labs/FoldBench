

from evaluation import eval_by_dockqv2,eval_by_ost
import pandas as pd
import argparse
import os



parser = argparse.ArgumentParser()
parser.add_argument(
    "--targets_dir", required=False, default='./examples/targets', help="The dir with the targets files."
)
parser.add_argument(
    "--evaluation_dir", required=False,default='./examples/outputs/evaluation', help="The dir with the evaluation files.",
)
parser.add_argument(
    "--algorithm_name", required=False, default='Protenix', help="The name of the algorithm.",
)
parser.add_argument(
    "--ground_truth_dir", required=False, default='./examples/ground_truths', help="The dir with the ground truth files.",
)

parser.add_argument(
        "--targets", required=False, default= ["interface_protein_ligand","interface_antibody_antigen","interface_protein_dna", "monomer_protein"], nargs='+', help="targets to evaluate.",
    )
args = parser.parse_args()

evaluation_dir = os.path.join(args.evaluation_dir,args.algorithm_name)

os.makedirs(os.path.join(evaluation_dir,'raw'), exist_ok=True)
target_types =  args.targets



prediction_summary_path = f'{evaluation_dir}/prediction_reference.csv'
prediction_summary_df = pd.read_csv(prediction_summary_path)

# caculation
for target_type in target_types:
    target_df_path = f'{args.targets_dir}/{target_type}.csv'
    if not os.path.exists(target_df_path):
        print(f"target_df_path is not exists for {target_type}")
        continue
    target_df = pd.read_csv(target_df_path)

    target_df = pd.merge(target_df,prediction_summary_df, on='pdb_id', how='left')

    if target_type in  ["interface_protein_protein","interface_antibody_antigen","interface_protein_peptide","interface_protein_ligand","interface_protein_dna","interface_protein_rna","monomer_dna","monomer_rna","monomer_protein"]:
        eval_by_ost(target_df,target_type,evaluation_dir,args.ground_truth_dir)
    
        if target_type in  ["interface_protein_dna","interface_protein_rna"]:
            eval_by_dockqv2(target_df,target_type,evaluation_dir,args.ground_truth_dir)