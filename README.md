# FoldBench: An All-atom Benchmark for Biomolecular Structure Prediction

![Abstract](./assets/fig1_abstract.png)


FoldBench is a low-homology benchmark that spans proteins, nucleic acids, ligands, and six major interaction types, enabling assessments that were previously infeasible with task-specific datasets.

## 🏆 Leaderboard

Scores represent the success rate for interface prediction tasks and LDDT for monomer prediction tasks.

### Protein Interactions

| Model | Protein-Protein | Antibody–Antigen | Protein-Ligand |
|:--------------:|:--------------:|:-----:|:--------------:|
| AlphaFold 3    | 72.92%          | 47.90%  | 64.90%           |
| Boltz-1        | 68.25%          | 33.54% | 55.04%          |
| Chai-1         | 68.53%          | 23.64% | 51.23%          |
| HelixFold 3    | 66.27%          | 28.40%  | 51.82%          |
| Protenix       | 68.18%          | 34.13% | 50.70%           |

### Nucleic acids

| Model | Protein-RNA | Protein-DNA | RNA Monomer | DNA Monomer |
|:--------------:|:-----------:|:-----------:|:-----------:|:-----------:|
| AlphaFold 3    | 62.32%       | 79.18%       | 61.40%        | 52.78%        |
| Boltz-1        | 56.90%       | 70.97%       | 44.16%        | 33.88%        |
| Chai-1         | 50.91%       | 69.97%       | 48.52%        | 45.74%        |
| HelixFold 3    | 48.28%       | 50.00%       | 54.87%        | 29.09%        |
| Protenix       | 44.78%       | 68.39%       | 59.03%        | 43.84%        |

**Note:**
- Interface prediction is evaluated by success rate.
- Monomer prediction is evaluated by LDDT.
- Success is defined as:
  - For protein–ligand interfaces: LRMSD < 2 Å and LDDT-PLI > 0.8
  - For all other interfaces: DockQ ≥ 0.23

## 🚀 Getting Started

To get started with FoldBench, clone the repository and set up the Conda environment.

```bash
# 1. Clone the repository
git clone https://github.com/BEAM-Labs/FoldBench.git
cd FoldBench

# 2. Create and activate the Conda environment for evaluation
conda env create -f environment.yml
conda activate foldbench
```

## ⚙️ Reproducing Our Evaluation

You can use our evaluation scripts to reproduce the scores for any model, such as Protenix, whose predictions are included as an example. The final results will be generated in `summary_table.csv`.

```bash
# Ensure you are in the FoldBench root directory and the conda environment is active

# Step 1: Calculate per-target scores from prediction files
# This uses OpenStructure (ost) and DockQ to score each prediction against its ground truth
python evaluate.py \
  --targets_dir ./targets \
  --evaluation_dir ./outputs/evaluation \
  --algorithm_name Protenix \
  --ground_truth_dir ./ground_truths

# Step 2: Aggregate scores and calculate the final success rates/LDDT
# This summarizes the results for specified models and tasks into a final table
python task_score_summary.py \
  --evaluation_dir ./outputs/evaluation \
  --target_dir ./targets \
  --output_path ./summary_table.csv \
  --algorithm_names Protenix \
  --targets interface_protein_ligand interface_protein_dna monomer_protein \
  --metric_type rank
```

## ✨ Integrating a New Model into FoldBench

We enthusiastically welcome community submissions\!

To ensure a fair and unbiased comparison, FoldBench operates as a **blind benchmark**. We do not publicize the target PDBs. Instead, you can submit your algorithm for us to run the tests.

For detailed instructions on how to package your model for submission, please see the contributor's guide:
**[Integrating a New Model into FoldBench](./algorithms/README.md)**

## 🙏 Acknowledgements

We gratefully acknowledge the developers of the following projects, which are essential to FoldBench:

+ [Alphafold3](https://github.com/google-deepmind/alphafold3)
+ [Protenix](https://github.com/bytedance/Protenix)
+ [Chai-1](https://github.com/chaidiscovery/chai-lab)
+ [Boltz-1](https://github.com/jwohlwend/boltz)
+ [Helixfold3](https://github.com/PaddlePaddle/PaddleHelix/tree/dev/apps/protein_folding/helixfold3)
+ [OpenStructure](https://git.scicore.unibas.ch/schwede/openstructure)
+ [DockQ](https://github.com/bjornwallner/DockQ)

## ✍️ How to Cite

If you use FoldBench in your research, please cite our paper:

```bibtex
@article {Xu2025.05.22.655600,
	author = {Xu, Sheng and Feng, Qiantai and Qiao, Lifeng and Wu, Hao and Shen, Tao and Cheng, Yu and Zheng, Shuangjia and Sun, Siqi},
	title = {FoldBench: An All-atom Benchmark for Biomolecular Structure Prediction},
	year = {2025},
	doi = {10.1101/2025.05.22.655600},
	journal = {bioRxiv}
}
```


