# TCGA Lung Methylation Analysis

This repository contains parallel DNA methylation analysis workflows for The
Cancer Genome Atlas lung adenocarcinoma (`TCGA-LUAD`) and lung squamous cell
carcinoma (`TCGA-LUSC`) cohorts.

The scripts download Illumina HumanMethylation450 data, prepare methylation and
clinical objects, identify differentially methylated CpGs (DMCs), and detect
and annotate differentially methylated regions (DMRs).

## Repository layout

```text
analysis/
├── luad/   # TCGA-LUAD methylation workflow
└── lusc/   # TCGA-LUSC methylation workflow
```

Both directories use the same numbered workflow:

1. `1DataDownload.R` queries and downloads TCGA methylation data.
2. `2Preparation.R` prepares and filters methylation matrices.
3. `3DMCs.R` tests individual CpG sites with `limma`.
4. `4DMRs.R` identifies and annotates DMRs with `DMRcate`.
5. `5DMRs.R` provides additional visualization and interpretation.

Run the scripts in numerical order for the cancer type being analyzed. These
are historical analysis scripts rather than a fully automated pipeline.

## Environment

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate tcga-lung-methylation
```

## Input locations

The scripts use environment variables instead of machine-specific paths:

- `LUAD_DATA_DIR`: working directory for LUAD inputs and generated files;
- `LUSC_DATA_DIR`: working directory for LUSC inputs and generated files;
- `CPG_ISLAND_BED`: CpG-island annotation BED file;
- `DNASE_SITES_BED`: DNase-site annotation BED file;
- `CROSS_REACTIVE_PROBES_FILE`: cross-reactive 450K probe annotation CSV.

For example:

```bash
export LUAD_DATA_DIR=/path/to/controlled-data/luad
export LUSC_DATA_DIR=/path/to/controlled-data/lusc
```

When an annotation variable is unset, the scripts look for the corresponding
file in the active data directory.

## Data policy

Downloaded GDC files, RDS objects, clinical tables, methylation matrices,
sample identifiers, intermediate objects, and generated results are excluded
from Git. Store any controlled or participant-level data in an approved data
location outside this repository.

The repository contains analysis code only. Review the applicable GDC/TCGA
data-access terms before downloading or redistributing data.
