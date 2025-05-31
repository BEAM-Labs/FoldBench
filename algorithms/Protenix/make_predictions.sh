#!/bin/bash

af3_input_json=$1
input_dir=$2
prediction_dir=$3
evaluation_dir=$4
gpu_id=$5

PYTHON_PATH="/opt/conda/bin/python"
# convert af3 input data to model format
$PYTHON_PATH ./preprocess.py --af3_input_json="$af3_input_json" --input_dir="$input_dir" 


# run inference

export CUDA_VISIBLE_DEVICES=$gpu_id
N_sample=5
N_step=200
N_cycle=10

# seed=42
seed=42,66,101,2024,8888

$PYTHON_PATH /algo/Protenix/runner/inference.py \
--seeds ${seed} \
--dump_dir ${prediction_dir} \
--input_json_path $input_dir/inputs.json \
--model.N_cycle ${N_cycle} \
--sample_diffusion.N_sample ${N_sample} \
--sample_diffusion.N_step ${N_step}  \
--use_msa_server

# Convert predictions to the general cif format, 
# and generate evaluation prediction_summary.csv in evaluation_dir
$PYTHON_PATH ./postprocess.py --input_dir="$input_dir" --prediction_dir="$prediction_dir" --evaluation_dir="$evaluation_dir"