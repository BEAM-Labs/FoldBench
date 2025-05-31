"""
preprocess for algorithm. 
For example, you generally need to convert AF3 input JSON files to the algorithm expected format.
You can use PreProcess.preprocess() to perform preprocessing.
"""

import argparse
import os
import json
from tqdm import tqdm


class PreProcess():
    def __init__(self):
        super().__init__()
    
    def format_single_input(self, input_data):
        tmp_seq_data={}
        tmp_seq_data['name']= input_data['name']
        tmp_seq_data["sequences"]=[]

        for seq in input_data['sequences']: # seq is a dict
            tmp_seq={}

            seq_type = next(iter(seq))

            if seq_type=='protein':
                tmp_seq["proteinChain"]={}
                tmp_seq["proteinChain"]['sequence']=seq[seq_type]['sequence']
                if isinstance(seq[seq_type]['id'],list):
                    tmp_seq["proteinChain"]['count']=len(seq[seq_type]['id'])
                else:
                    tmp_seq["proteinChain"]['count']=1
                
                
                # add modification
                tmp_seq["proteinChain"]['modifications']=[]
                for modification in seq[seq_type]['modifications']:
                    tmp_modification = {}
                    tmp_modification["ptmType"]=f'CCD_{modification["ptmType"]}'
                    tmp_modification['ptmPosition']=modification['ptmPosition']
                    tmp_seq["proteinChain"]['modifications'].append(tmp_modification)
                

            elif seq_type=='rna':
                tmp_seq["rnaSequence"]={}
                tmp_seq["rnaSequence"]['sequence']=seq[seq_type]['sequence']
                if isinstance(seq[seq_type]['id'],list):
                    tmp_seq["rnaSequence"]['count']=len(seq[seq_type]['id'])
                else:
                    tmp_seq["rnaSequence"]['count']=1
                
                # add modification
                tmp_seq["rnaSequence"]['modifications']=[]
                for modification in seq[seq_type]['modifications']:
                    tmp_modification = {}
                    tmp_modification["modificationType"]=f'CCD_{modification["modificationType"]}'
                    tmp_modification['basePosition']=modification['basePosition']
                    tmp_seq["rnaSequence"]['modifications'].append(tmp_modification)

            elif seq_type=='dna':
                tmp_seq["dnaSequence"]={}
                tmp_seq["dnaSequence"]['sequence']=seq[seq_type]['sequence']
                if isinstance(seq[seq_type]['id'],list):
                    tmp_seq["dnaSequence"]['count']=len(seq[seq_type]['id'])
                else:
                    tmp_seq["dnaSequence"]['count']=1
                
                # add modification
                tmp_seq["dnaSequence"]['modifications']=[]
                for modification in seq[seq_type]['modifications']:
                    tmp_modification = {}
                    tmp_modification["modificationType"]=f'CCD_{modification["modificationType"]}'
                    tmp_modification['basePosition']=modification['basePosition']
                    tmp_seq["dnaSequence"]['modifications'].append(tmp_modification)
                
            elif seq_type=='ligand':
                tmp_seq["ligand"]={}
                
                tmp_seq["ligand"]['ligand']="CCD"
                for ccd in seq[seq_type]['ccdCodes']:
                    tmp_seq["ligand"]['ligand'] += f'_{ccd}'
                
                if isinstance(seq[seq_type]['id'],list):
                    tmp_seq["ligand"]['count']=len(seq[seq_type]['id'])
                else:
                    tmp_seq["ligand"]['count']=1

            tmp_seq_data["sequences"].append(tmp_seq)

        return tmp_seq_data

    def preprocess(self, af3_input_json_path, input_dir):
        """
        af3_input_json_path: path to the af3 input json file
        input_dir: path to the input directory for the algorithm
        """

        with open(af3_input_json_path, "r") as f:
            folding_inputs = json.load(f)
        
        
        mapped_folding_inputs = [
            self.format_single_input(folding_inputs[i])
            for i in tqdm(range(len(folding_inputs)))
        ]
    
        with open(f'{input_dir}/inputs.json', "w") as f:
            json.dump(mapped_folding_inputs, f,indent = 4)
        print(
            "{} folding inputs written to {}.".format(len(mapped_folding_inputs), input_dir)
        )
           

parser = argparse.ArgumentParser()
parser.add_argument(
    "--af3_input_json",
    help="The path to the input .json file.",
)
parser.add_argument(
    "--input_dir",
    help="The path to write prepared input data in the format expected by the algorithm.",
)
args = parser.parse_args()

preprocess = PreProcess()
preprocess.preprocess(args.af3_input_json,args.input_dir)
