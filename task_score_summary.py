import pandas as pd
import os
import glob
import argparse

metric_max = ['dockq_score','lddt-lp','lddt-pli','gdt-ts','tm-score','lddt']
metric_min = ['irmsd','lrmsd','rmsd']

target_metrics_summary = {}
target_metrics_summary['interface_protein_ligand'] = {'rmsd_lddt-pli_success_rate': [], 'lddt-lp': [], 'lddt-pli': []}
target_metrics_summary['interface_protein_protein'] = {'dockq_score_success_rate': [], 'irmsd': [], 'lrmsd': [], 'lddt': []}
target_metrics_summary['interface_antibody_antigen'] = {'dockq_score_success_rate': [], 'irmsd': [], 'lrmsd': [], 'lddt': []}
target_metrics_summary['interface_protein_dna'] = {'dockq_score_success_rate': [], 'irmsd': [], 'lrmsd': [], 'lddt': []}
target_metrics_summary['interface_protein_rna'] = {'dockq_score_success_rate': [], 'irmsd': [], 'lrmsd': [], 'lddt': []}
target_metrics_summary['interface_protein_peptide'] = {'dockq_score_success_rate': [], 'irmsd': [], 'lrmsd': [], 'lddt': []}
target_metrics_summary['monomer_protein'] = {'gdt-ts': [], 'tm-score': [], 'rmsd': [], 'lddt': []}
target_metrics_summary['monomer_dna'] = {'gdt-ts': [], 'tm-score': [], 'rmsd': [], 'lddt': []}
target_metrics_summary['monomer_rna'] = {'gdt-ts': [], 'tm-score': [], 'rmsd': [], 'lddt': []}


def change_column_name(df):
    if 'native_chain_id_1' in df.columns and 'native_chain_id_2' in df.columns:
        df.rename(columns={'native_chain_id_1': 'interface_chain_id_1', 'native_chain_id_2': 'interface_chain_id_2'}, inplace=True)
    return df

def find_overlap_sample(df, ref_df):

    if 'interface_chain_id_1' in df.columns and 'interface_chain_id_2' in df.columns:
        df_tuples = list(zip(df['pdb_id'], 
                        df['interface_chain_id_1'], 
                        df['interface_chain_id_2']))


        ref_df_tuples = list(zip(ref_df['pdb_id'], 
                                ref_df['interface_chain_id_1'], 
                                ref_df['interface_chain_id_2']))
    else:
        df_tuples = list(zip(df['pdb_id']))
        ref_df_tuples = list(zip(ref_df['pdb_id']))

    

    common_tuples = set(df_tuples).intersection(set(ref_df_tuples))


    if 'interface_chain_id_1' in df.columns and 'interface_chain_id_2' in df.columns:
        mask = [tuple(row) in common_tuples for row in df[['pdb_id', 'interface_chain_id_1', 'interface_chain_id_2']].values]
    else:
        mask = [tuple(row) in common_tuples for row in df[['pdb_id']].values]

    
    df = df[mask]
    return df

# Get best prediction for each target
def get_best_rows(df, metric, metric_type):
    if metric_type == 'rank':
        if 'interface_chain_id_1' in df.columns and 'interface_chain_id_2' in df.columns:
            best_rows = df.loc[df.groupby(["pdb_id","interface_chain_id_1","interface_chain_id_2"])['ranking_score'].idxmax()]
        else:
            best_rows = df.loc[df.groupby(["pdb_id"])['ranking_score'].idxmax()]
    elif metric_type == 'best':
        if metric in metric_max:
            if 'interface_chain_id_1' in df.columns and 'interface_chain_id_2' in df.columns:
                best_rows = df.loc[df.groupby(["pdb_id","interface_chain_id_1","interface_chain_id_2"])[metric].idxmax()]
            else:
                best_rows = df.loc[df.groupby(["pdb_id"])[metric].idxmax()]
        elif metric in metric_min:
            if 'interface_chain_id_1' in df.columns and 'interface_chain_id_2' in df.columns:
                best_rows = df.loc[df.groupby(["pdb_id","interface_chain_id_1","interface_chain_id_2"])[metric].idxmin()]
            else:
                best_rows = df.loc[df.groupby(["pdb_id"])[metric].idxmin()]
    return best_rows

def calculate_success_rate(df, metric, metric_type):
    """
    Calculate success rate based on a metric and threshold.
    
    Args:
        df (pd.DataFrame): Input dataframe
        threshold (float): Threshold value for success
    Returns:
        float: Success rate
    """



    if metric  == 'dockq_score':
        df=df[df[metric].notna()]
        best_rows = get_best_rows(df, metric,metric_type)

        success_rate = (len(best_rows[best_rows[metric] >= 0.23])/len(best_rows))*100
        
    elif metric == 'rmsd':
        df=df[df[metric].notna()]
        best_rows = get_best_rows(df, metric,metric_type)
        success_rate = (len(best_rows[best_rows[metric] < 2.0])/len(best_rows))*100
    elif metric == 'rmsd_lddt-pli':
        df=df[df['rmsd'].notna() & df['lddt-pli'].notna()]
        best_rows = get_best_rows(df, metric,metric_type)
        success_rate = (len(best_rows[(best_rows['rmsd']<2.0) & (best_rows['lddt-pli']>0.8)])/len(best_rows))*100

    return success_rate


def calculate_score_avg(df, metric, metric_type):
    """
    Calculate success rate based on a metric and threshold.
    
    Args:
        df (pd.DataFrame): Input dataframe
        threshold (float): Threshold value for success
    Returns:
        float: Success rate
    """

    df=df[df[metric].notna()]
    best_rows = get_best_rows(df, metric,metric_type)
    avg_score = best_rows[metric].mean()


        
    return avg_score

def process_csv_files(evaluation_dir,target_dir,output_path,models,targets,metric_type):

    results = {}
    
    # Process interface files
    for model in models:
        results[model] = {}

        for target in targets:
            result_path = os.path.join(evaluation_dir,model,'raw',f"{target}_ost.csv")
            # 
            if not os.path.exists(result_path):
                print(f'{result_path} not found')
                continue
            target_path = os.path.join(target_dir, f"{target}.csv")

            result_df = pd.read_csv(result_path)
            result_df = change_column_name(result_df)

            target_df = pd.read_csv(target_path)
            target_df = change_column_name(target_df)
            print(f'{target} num: {len(target_df)}')

            result_df = find_overlap_sample(result_df,target_df)

            results[model][target] = {}

            # Process pp interface: dockq success_rate,irmsd,lrmsd,lddt
            if target in ["interface_protein_protein", "interface_protein_peptide", "interface_antibody_antigen","interface_protein_dna", "interface_protein_rna"]:
                   
                results[model][target]['lddt'] = calculate_score_avg(result_df, 'lddt',metric_type)
                
                
                if target in ["interface_protein_dna", "interface_protein_rna"]:
                    if os.path.exists(os.path.join(evaluation_dir,model,'raw', f"{target}_dockqv2.csv")):
                        result_path_dockqv2 = os.path.join(evaluation_dir,model,'raw', f"{target}_dockqv2.csv")
                        result_df_dockqv2 = pd.read_csv(result_path_dockqv2)
                        result_df_dockqv2 = change_column_name(result_df_dockqv2)
                        result_df_dockqv2 = find_overlap_sample(result_df_dockqv2,target_df)
                        results[model][target]['dockq_score_success_rate']  = calculate_success_rate(result_df_dockqv2, 'dockq_score',metric_type)
                        results[model][target]['irmsd'] = calculate_score_avg(result_df_dockqv2, 'irmsd',metric_type)
                        results[model][target]['lrmsd'] = calculate_score_avg(result_df_dockqv2, 'lrmsd',metric_type)
                else:
                    results[model][target]['dockq_score_success_rate']  = calculate_success_rate(result_df, 'dockq_score',metric_type)
                    results[model][target]['irmsd'] = calculate_score_avg(result_df, 'irmsd',metric_type)
                    results[model][target]['lrmsd'] = calculate_score_avg(result_df, 'lrmsd',metric_type)

            # Process pl interface: rmsd_lddt-pli success_rate,lddt-lp,lddt-pli
            elif target in ["interface_protein_ligand"]:
                results[model][target]['rmsd_lddt-pli_success_rate'] = calculate_success_rate(result_df, 'rmsd_lddt-pli',metric_type)
                results[model][target]['lddt-lp'] = calculate_score_avg(result_df, 'lddt-lp',metric_type)
                results[model][target]['lddt-pli'] = calculate_score_avg(result_df, 'lddt-pli',metric_type)

            # Process monomer: gdt_ts,tm-score,rmsd,lddt
            elif target in ["monomer_dna", "monomer_rna", "monomer_protein"]:
                if 'gdt_ts' in result_df.columns:
                    results[model][target]['gdt-ts'] = calculate_score_avg(result_df, 'gdt_ts',metric_type)
                elif 'gdt-ts' in result_df.columns:
                    results[model][target]['gdt-ts'] = calculate_score_avg(result_df, 'gdt-ts',metric_type)
                if 'tm-score' in result_df.columns:
                    results[model][target]['tm-score'] = calculate_score_avg(result_df, 'tm-score',metric_type)
                elif 'tm_score' in result_df.columns:
                    results[model][target]['tm-score'] = calculate_score_avg(result_df, 'tm_score',metric_type)
                
                results[model][target]['rmsd'] = calculate_score_avg(result_df, 'rmsd',metric_type)
                results[model][target]['lddt'] = calculate_score_avg(result_df, 'lddt',metric_type)
        
    return results

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--evaluation_dir", required=False, default='./outputs/evaluation', help="The dir with the evaluation files.",
    )
    parser.add_argument(
        "--target_dir", required=False, default='./targets', help="The dir with the target files.",
    )
    parser.add_argument(
        "--output_path", required=False, default='./summary_table.csv', help="output path",
    )
    parser.add_argument(
        "--algorithm_names", required=False, default= ['Protenix'], nargs='+', help="models to evaluate",
    )
    parser.add_argument(
        "--targets", required=False, default= ["interface_protein_ligand","interface_antibody_antigen","interface_protein_dna", "monomer_protein"], nargs='+', help="targets to evaluate.",
    )
    parser.add_argument(
        "--metric_type", required=False, default= "rank", help="rank or best",
    )
    args = parser.parse_args()
    results = process_csv_files(args.evaluation_dir,args.target_dir,args.output_path,args.algorithm_names,args.targets,args.metric_type)
    
    for target in args.targets:
        for metric in target_metrics_summary[target].keys():
            for model in args.algorithm_names:

                if target not in results[model].keys():
                    print(f'{model}: {target} not found')
                    continue
                if metric in results[model][target].keys():
                    target_metrics_summary[target][metric].append(results[model][target][metric])
                else:
                    target_metrics_summary[target][metric].append(None)

    # Create a DataFrame with models as columns
    df_wide = pd.DataFrame()
    for target in target_metrics_summary.keys():
        for metric in target_metrics_summary[target].keys():
            row_data = {
                'target': target,
                'metric': metric
            }
            # Add model values as columns
            for i, model in enumerate(args.algorithm_names):
                if i < len(target_metrics_summary[target][metric]) and target_metrics_summary[target][metric][i] is not None:
                    row_data[model] = round(target_metrics_summary[target][metric][i],2)
                else:
                    row_data[model] = None
            if row_data[model] is not None:
                df_wide = pd.concat([df_wide, pd.DataFrame([row_data])], ignore_index=True)
    
    # Save the wide format DataFrame
    df_wide.to_csv(args.output_path, index=False)

    # Print the table in a nicely formatted way
    print("\nResults Summary:")
    print("=" * 80)
    
    # Group by target and print each group
    for target in df_wide['target'].unique():
        print(f"\n{target.upper()}")
        print("-" * 80)
        
        # Get rows for this target
        target_df = df_wide[df_wide['target'] == target]
        
        # Print header
        header = "Metric".ljust(20)
        for model in args.algorithm_names:
            header += f"{model}".rjust(15)
        print(header)
        print("-" * 80)
        
        # Print each metric row
        for _, row in target_df.iterrows():
            metric_str = str(row['metric']).ljust(20)
            for model in args.algorithm_names:
                value = row[model]
                if pd.isna(value):
                    metric_str += "N/A".rjust(15)
                else:
                    metric_str += f"{value:.2f}".rjust(15)
            print(metric_str)
    
    print("\n" + "=" * 80)
