# Nanopore Basecalling Benchmarking Pipeline (AUC)

Pipeline for benchmarking Oxford Nanopore basecalling-calling performance across different basecallers, models, and q-score thresholds using ROC AUC analysis. 

## Repository Structure

```bash id="bov8os"
repo/
├── run_dorado_fast_new_EpiC.sh
├── run_dorado_hac_new_EpiC.sh
├── run_dorado_hac_new_EpiC_server.sh
├── run_dorado_sup_new_EpiC.sh
├── run_guppy_hac_new_EpiC.sh
├── align_methyl_bams.sh
├── auc_model_qscore_analysis.py
├── plot_auc_results.py
└── README.md
```

## Scripts

| Script                              | Purpose                                         |
| ----------------------------------- | ----------------------------------------------- |
| `run_dorado_fast_new_EpiC.sh`       | Dorado FAST methylation calling                 |
| `run_dorado_hac_new_EpiC.sh`        | Dorado HAC methylation calling                  |
| `run_dorado_hac_new_EpiC_server.sh` | Dorado HAC methylation calling (server version) |
| `run_dorado_sup_new_EpiC.sh`        | Dorado SUP methylation calling                  |
| `run_guppy_hac_new_EpiC.sh`         | Guppy HAC methylation calling                   |
| `align_methyl_bams.sh`              | Align, sort, and index BAM files                |
| `auc_model_qscore_analysis.py`      | Compute AUC benchmarking metrics                |
| `plot_auc_results.py`               | Generate benchmarking plots                     |

## Input Directory Structure

project/
│
├── pod5_data/
│   ├── sample_1/
│   ├── sample_2/
│   └── ...
│
├── fast5_data/
│   ├── sample_1/
│   ├── sample_2/
│   └── ...
│
├── reference.fasta
├── reference.fasta.fai
│
├── run_dorado_fast_new_EpiC.sh
├── run_dorado_hac_new_EpiC.sh
├── run_dorado_hac_new_EpiC_server.sh
├── run_dorado_sup_new_EpiC.sh
├── run_guppy_hac_new_EpiC.sh
│
├── align_methyl_bams.sh
├── auc_model_qscore_analysis.py
├── plot_auc_results.py
│


## Workflow

### 1. Generate Methylation BAMs

```bash id="i4h25f"
bash run_dorado_fast_new_EpiC.sh
bash run_dorado_hac_new_EpiC.sh
bash run_dorado_sup_new_EpiC.sh
bash run_guppy_hac_new_EpiC.sh
```

### 2. Align and Sort BAMs

```bash id="4m6g7x"
bash align_methyl_bams.sh <bam_root_dir> <reference_fasta> <output_dir> 8
```

Example:

```bash id="m9rx1q"
bash align_methyl_bams.sh aligned_methyl reference.fasta aligned_methyl_7Apr 8
```

### 3. Run AUC Benchmarking

```bash id="lnh6eq"
python auc_model_qscore_analysis.py aligned_methyl_7Apr reference.fasta
```

### 4. Generate Plots

```bash id="t22n7d"
python plot_auc_results.py auc_results
```

## Main Outputs

* BAM files
* Sorted/indexed BAMs
* BED cache files
* `model_qscore_summary.csv`
* `per_sample_auc.csv`
* `top5_combinations.csv`
* AUC benchmarking plots

## Requirements

* Python 3.9+
* Dorado
* Guppy
* modkit
* minimap2
* samtools

## Notes

* Edit input/output/reference paths inside scripts before running.

## Full Documentation

[Google Docs Documentation](https://docs.google.com/document/d/1AxCqu6Y9VaVnnh5uNfBXIIr3nSeZ_k0oSflzP23R1LA/edit?tab=t.0)
