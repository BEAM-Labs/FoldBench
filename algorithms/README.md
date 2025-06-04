# Integrating a New Model into FoldBench

Welcome, contributor! This guide outlines the process for integrating your protein structure prediction model into the FoldBench benchmarking platform.

## Platform Requirements

Before you begin, ensure your environment meets the following system requirements:
* **OS:** Linux
* **Containerization:** Apptainer
* **Package Management:** Conda
* **Hardware:** NVIDIA GPU

---

## Preparing Your Algorithm

To add your model, create a new directory inside `FoldBench/algorithms/` using your algorithm's name. This directory must contain the following four files:

### 1. 📦 `container.def`
This Apptainer definition file specifies the complete environment for your model. It should install all necessary system packages, Python libraries (e.g., via pip), and set up any required environment variables.

### 2. 📝 `preprocess.py`
This script prepares the input data for your model. It must contain a `PreProcess.preprocess()` method.

* **Input:** The script should read from the benchmark's standard input file: `./alphafold3_inputs.json`.
* **Function:** Convert the data from the standard JSON format into the specific file format that your model requires for inference.
* **Output:** The preprocessed data should be saved to your algorithm's dedicated input directory: `./outputs/input/{algorithm_name}/`.

### 3. 🚀 `make_predictions.sh`
This is the main inference script that runs your model. It will be executed from within the Apptainer environment.

* **Input:** It should read the preprocessed data from `./outputs/input/{algorithm_name}/`.
* **Function:** Execute your model's prediction command-line interface.
* **Output:** The prediction artifacts (e.g., `.cif` or `.pdb` files) must be written to the prediction directory: `./outputs/prediction/{algorithm_name}/`.

### 4. ✨ `postprocess.py`
This script standardizes your model's output for evaluation. It must contain a `PostProcess.postprocess()` method and perform two key tasks:

1.  **Generate Prediction Summary:** Create a summary file named `prediction_reference.csv` in the evaluation directory: `./outputs/evaluation/{algorithm_name}/prediction_reference.csv`. This CSV file is **required** for the benchmark and must include the following columns: `pdb_id`, `seed`, `sample`, `ranking_score`, and `prediction_path`.
2.  **Format for Evaluation:** Convert your model's raw output files (located in `./outputs/prediction/{algorithm_name}/`) into a format compatible with our evaluation tools ([OpenStructure](https://git.scicore.unibas.ch/schwede/openstructure) and [DockQ](https://github.com/bjornwallner/DockQ)).

---

## Running the Benchmark

Once you have prepared your four files, you can test the entire workflow using our provided scripts.

### Step 1: Build Environments
Navigate to the root `FoldBench/` directory and run the build script. This command builds the Apptainer image for your algorithm. 

```bash
cd FoldBench

# Build the Apptainer image for your algorithm
./build_apptainer_images.sh

# Create the conda environment for evaluation
conda env create -f environment.yml
```

### Step 2: Run Inference and Evaluation
Activate the conda environment and execute the main run script. This will automate the preprocessing, prediction, postprocessing, and scoring for all registered algorithms.

```bash
conda activate foldbench
./run.sh
```

Once your model runs successfully, please submit a pull request to add it to our platform. We look forward to your contribution!