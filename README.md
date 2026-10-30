# Hurst exponent estimation and analysis

Code for **Hurst exponent estimation and analysis** of fMRI and eye-tracking signals, including associations with MRS/fMRS measurements, cortical cytoarchitecture, and functional gradients. The repository contains Python, MATLAB, and R scripts, Jupyter notebooks, and code used to produce the manuscript figures.

## Data availability

Raw data are not included in this repository:

- **Discovery (WAND):** [GIN dataset and description](https://doi.org/10.1038/s41597-024-04154-7).
- **Replication:** [Processed, anonymized data on Zenodo](https://doi.org/10.5281/zenodo.15855362); additional information is available from [Prof. Jorge Jovicich](mailto:jorge.jovicich@unitn.it) upon reasonable request.
- **Robustness (HCP Young Adult):** [100-unrelated sample, ConnectomeDB](https://db.humanconnectome.org).

Preprocessing and Hurst estimation scripts require locally configured input/output paths and external data. Downstream analyses use the derived files in `DATA/` and `RESULTS/`.

## Repository structure

```text
CODE/
├── datasets/
│   ├── discovery_dataset/
│   ├── replication_dataset/
│   └── robustness_dataset/
├── eye_tracking/
├── figures/
│   ├── main/
│   └── supplementary/
└── utils/
DATA/
├── cytoarchitecture/
├── discovery_dataset/
├── gradients/
├── pupillometry/
├── replication_dataset/
├── robustness_dataset/
└── templates/
RESULTS/
├── discovery_dataset/
├── replication_dataset/
├── robustness_dataset/
├── pupillometry/
└── spin_tests_correlations/
FIGURES/
├── main/
└── supplementary/
```

### `CODE/`

- **`datasets/discovery_dataset/`** — WAND resting-state fMRI and MRS: parcellation (`01_preprocessing/`), Hurst estimation (`02_hurst_estimation/`), MRS quality control and regional statistics (`03_analyses/`), and models relating Hurst to spectroscopy, cytoarchitecture, and gradients (`04_models/`). `alternative_parcellation_glasser360/` repeats selected steps with Glasser360.
- **`datasets/replication_dataset/`** — Resting-state and working-memory fMRI analyses: parcellation, Hurst estimation, task/network statistics and SVD, MRS/fMRS models, functional-gradient models, and Glasser360 analyses. Includes an FSL script for voxelwise group inference.
- **`datasets/robustness_dataset/`** — HCP resting-state Hurst summaries and functional-gradient models, using previously estimated Hurst values.
- **`eye_tracking/`** — Pupil-size and microsaccade preprocessing, removal of task-evoked components, Hurst estimation, dataframe preparation, and statistical models of session order and neurometabolite measures.
- **`figures/main/` and `figures/supplementary/`** — Notebooks generating manuscript figures from saved data and model results.
- **`utils/`** — Shared functions for parcellation, statistics, cortical surface plots, network profiles, SVD, and contrasts.

Within the dataset and eye-tracking pipelines, numbered subfolders follow the usual sequence: `01_preprocessing/` → `02_hurst_estimation/` → `03_analyses/` → `04_models/` (where applicable).

### `DATA/`, `RESULTS/`, and `FIGURES/`

- **`DATA/discovery_dataset/`** — WAND MRS tables, quality-control information, and MRS voxel-to-parcel correspondence.
- **`DATA/replication_dataset/`** — MRS/fMRS-derived metabolite tables and voxel-to-parcel correspondence.
- **`DATA/robustness_dataset/`** — HCP-derived Hurst estimates.
- **`DATA/pupillometry/`** — Eye-tracking processing inputs and derived signals.
- **`DATA/cytoarchitecture/`** — von Economo classes and BigBrain microstructural gradient.
- **`DATA/gradients/`** — Cortical functional-gradient data.
- **`DATA/templates/`** — Cortical atlases, parcel labels, centroids, and surface annotations required by some analyses; obtain missing templates separately.
- **`RESULTS/`** — Hurst outputs, analysis-ready tables, statistical tests, mixed-model outputs, SVD results, and spatial permutation results, arranged by dataset or modality.
- **`FIGURES/main/` and `FIGURES/supplementary/`** — Exported figure files; some are omitted to reduce repository size.

## Requirements

- **Python 3.11** and **R 4.4** with packages specified in `environment.yml`.
- **MATLAB R2022b (9.13)** for Hurst estimation, with the [nonfractal](https://github.com/wonsang/nonfractal) and [WMTSA](https://www.atmos.washington.edu/~wmtsa/) toolboxes.
- **FSL** for the optional voxelwise permutation analysis; compatible EyeLink-reading dependencies for eye-tracking preprocessing.

Install the Conda environment from the repository root:

```bash
conda env create -f environment.yml
conda activate hurst-neuroimaging
```

## Running the analyses

1. **To preprocess data and estimate Hurst:** provide the appropriate fMRI/eye-tracking inputs and atlas files; configure paths in `01_preprocessing/` and `02_hurst_estimation/`, including MATLAB toolbox paths. The fMRI scripts parcellate *previously preprocessed* time series; raw scanner-level fMRI preprocessing is not included.
2. **To reproduce downstream analyses:** use the available intermediate files under `DATA/` and `RESULTS/` and run the relevant `03_analyses/` notebooks, `04_models/` R scripts, and `CODE/figures/` notebooks. Run R scripts from the repository root. Missing intermediate files or atlas templates must be supplied before running analyses that depend on them.

The fMRI Hurst scripts use wavelet-based maximum-likelihood estimation (`bfn_mfin_ml`, D4 MODWT, scales 1–5). Check sampling and estimation settings before applying the pipeline to another dataset.