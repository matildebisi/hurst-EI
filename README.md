# Hurst exponent analysis of fMRI, MRS, and eye-tracking data

Code and analysis outputs for a multimodal  study of long-range temporal dynamics, cortical organization, and their relationship to neurometabolite measurements. The repository combines regional fMRI Hurst exponent estimation with magnetic resonance spectroscopy (MRS), functional MRS (fMRS)-derived measurements, pupillometry, microsaccades, cortical cytoarchitecture, and macroscale functional gradients. It includes discovery, replication, and independent robustness analyses; an alternative cortical parcellation; and scripts/notebooks used to generate main and supplementary figures.

> **Data availability and scope.** Raw participant-level fMRI, MRS/fMRS, and eye-tracking recordings are **not distributed** with this repository. The preprocessing and Hurst estimation code is provided to document the complete methodological workflow and to allow use with appropriately obtained data. Those scripts require users to supply their own input files, atlas resources, external tools, and machine-specific paths. The downstream analyses are organized around the intermediate CSV/MAT files and saved results in `DATA/` and `RESULTS/`, allowing the main statistical and visualization workflow to be examined and reproduced without the original raw acquisitions **when all referenced intermediate inputs and templates are available**. Some datasets, template files, result files, and rendered figures may be excluded from a particular archive or GitHub release because of size or sharing restrictions; see [Included and missing files](#included-and-missing-files).

## Contents

1. [Scientific overview](#scientific-overview)
2. [Analysis workflow](#analysis-workflow)
3. [Repository structure](#repository-structure)
4. [Code: dataset-specific pipelines](#code-dataset-specific-pipelines)
5. [Code: eye tracking](#code-eye-tracking)
6. [Code: figures and shared utilities](#code-figures-and-shared-utilities)
7. [Data, results, and figures](#data-results-and-figures)
8. [Installation and dependencies](#installation-and-dependencies)
9. [How to reproduce the analyses](#how-to-reproduce-the-analyses)
10. [Data formats and naming conventions](#data-formats-and-naming-conventions)
11. [Included and missing files](#included-and-missing-files)
12. [Reproducibility, limitations, and responsible data sharing](#reproducibility-limitations-and-responsible-data-sharing)
13. [Citation and license](#citation-and-license)

## Scientific overview

### Main scientific questions

The analyses investigate the scale-free temporal structure of neurophysiological signals, quantified by the **Hurst exponent (H)**, and assess its dependence on cognitive state and cortical organization. The main components are:

- **Cortical fMRI dynamics:** regional Hurst exponents calculated from preprocessed parcelwise fMRI residual time series, summarized over the cortex, functional networks, participants, and experimental conditions.
- **Cognitive-load effects:** resting-state versus 0-back, 1-back, and 2-back conditions in the task-based dataset, with parcelwise and network-level statistical comparisons.
- **MRS and fMRS:** relationships of fMRI Hurst exponents to spectroscopy-derived GABA, Glx, and an excitation/inhibition proxy (`EIB`, based on Glx/GABA). The code includes spectroscopy quality control, task/run matching, subject-level analyses, and mixed-effects modeling. The replication input table includes measurements associated with static MRS and fMRS acquisitions; spectroscopy acquisition and spectral fitting themselves are not implemented in this repository.
- **Cortical architecture:** associations with von Economo cytoarchitectural classes and a BigBrain-derived microstructural gradient.
- **Functional organization:** relationships with the principal Margulies functional gradient, including polynomial trends and task-dependent interactions.
- **Peripheral physiological signals:** preprocessing and Hurst estimation of pupil-size and microsaccade time series, followed by models of task/session order and spectroscopy measures.
- **Robustness and sensitivity:** independent resting-state dataset, Schaefer100 versus Glasser360 parcellations, dimensionality-reduction analyses, and spatially informed statistical tests.

### Datasets and modalities

| Component | Role | Key modalities | Primary analysis |
|---|---|---|---|
| `discovery_dataset` | Initial analyses, including the WAND cohort | Resting-state fMRI; MRS | Hurst estimation, MRS quality control, spectroscopy-voxel comparisons, cytoarchitecture and gradients |
| `replication_dataset` | Task-based replication analyses | Rest/task fMRI; MRS and fMRS-derived metabolite table | Task and cortical network differences; GABA, Glx, and EIB models; gradient models; spatial statistics |
| `robustness_dataset` | Independent robustness assessment (HCP resting-state data) | Resting-state fMRI, provided here as derived estimates | Regional Hurst summaries and principal-gradient relationships |
| `eye_tracking` | Peripheral physiological analyses | Pupil diameter and microsaccades | Signal preprocessing, task-effect residualization, Hurst estimation, task-order and spectroscopy associations |
| `alternative_parcellation_glasser360` | Parcellation sensitivity analyses | fMRI | Repeat selected analyses using Glasser360 instead of Schaefer100 |

**Terminology:** `Schaefer100` is the 100-parcel cortical atlas; `Glasser360` is the 360-region HCP-MMP parcellation. `MRS` and `fMRS` refer to magnetic resonance spectroscopy and functional MRS, respectively. `Glx` represents the combined glutamate/glutamine signal, `GABA` the GABA measurement, and `EIB` the study's spectroscopy-derived Glx/GABA ratio or corresponding E/I proxy. This ratio is an *indirect proxy*, not a direct measurement of synaptic excitation or inhibition.

## Analysis workflow

The numbered subfolders indicate the **scientific sequence**, not a single executable build command:

```text
RAW DATA (not distributed)
    |
    +-- fMRI --> preprocessed voxelwise time series
    |              --> cortical parcellation [01_preprocessing]
    |              --> parcelwise Hurst estimation [02_hurst_estimation]
    |              --> per-participant .mat results
    |              --> analysis-ready CSV tables [03_analyses]
    |              --> statistics and mixed-effects models [03_analyses / 04_models]
    |              --> figures [CODE/figures]
    |
    +-- EyeLink eye tracking --> pupil / microsaccade preprocessing
    |                           --> task-evoked component residualization
    |                           --> Hurst estimation
    |                           --> dataframes and models --> figures
    |
    +-- MRS / fMRS --> spectroscopy processing outside this repository
                        --> provided metabolite tables and MRS QC
                        --> voxel/parcel matching and models --> figures

ADDITIONAL REFERENCE DATA:
    cortical atlases, MRS voxel-to-parcel correspondence,
    cytoarchitecture classes, BigBrain gradient, Margulies gradients
```

**Key separation:** stages `01_preprocessing/` and `02_hurst_estimation/` are supplied as methodological code and **must be adapted** to an individual site's file layout. The downstream notebooks and R scripts primarily use paths relative to the repository root. Their inputs are the processed tables and outputs represented by `DATA/` and `RESULTS/`; they do not require re-running raw preprocessing when the necessary derived files have been supplied.

### Hurst estimation settings

The MATLAB fMRI Hurst estimators call `bfn_mfin_ml` from the external `nonfractal` toolbox to obtain maximum-likelihood estimates. Their documented analysis configuration includes:

| Setting | Value used in the scripts |
|---|---|
| Estimation | Maximum likelihood (`ML`) via `bfn_mfin_ml` |
| Wavelet | Maximal overlap discrete wavelet transform (`modwt`) |
| Filter | `D4` |
| Scale range | `[1, 5]` |
| Detrending | Linear detrending and mean centering for fMRI |
| Additional decimation | Disabled in the main fMRI scripts (`decimate_b = false`) |
| Padding | Disabled in the main fMRI scripts (`padding = false`) |
| Output | Parcelwise Hurst estimates in MATLAB `.mat` files |

These are the parameters implemented in the scripts, not instructions to apply the settings blindly to a different acquisition protocol. **Pupillometry uses a separate Hurst estimation script** and input time-series preparation (including the stored `decimated_2s` configuration). Review the exact toolbox calls and relevant acquisition parameters before using another dataset.

## Repository structure

```text
hurst-analysis/
├── README.md                     # This document
├── environment.yml               # Conda environment with Python and R dependencies
├── .gitignore                    # Git exclusions (review before public release)
├── CODE/
│   ├── datasets/
│   │   ├── discovery_dataset/
│   │   │   ├── 01_preprocessing/
│   │   │   ├── 02_hurst_estimation/
│   │   │   ├── 03_analyses/
│   │   │   ├── 04_models/
│   │   │   └── alternative_parcellation_glasser360/
│   │   ├── replication_dataset/
│   │   │   ├── 01_preprocessing/
│   │   │   ├── 02_hurst_estimation/
│   │   │   ├── 03_analyses/
│   │   │   ├── 04_models/
│   │   │   └── alternative_parcellation_glasser360/
│   │   └── robustness_dataset/
│   │       ├── 03_analyses/
│   │       └── 04_models/
│   ├── eye_tracking/
│   │   ├── 01_preprocessing/
│   │   ├── 02_hurst_estimation/
│   │   ├── 03_analyses/
│   │   └── 04_models/
│   ├── figures/
│   │   ├── main/
│   │   └── supplementary/
│   └── utils/                   # Reusable Python functions
├── DATA/
│   ├── cytoarchitecture/
│   ├── discovery_dataset/
│   ├── gradients/
│   ├── pupillometry/
│   ├── replication_dataset/
│   ├── robustness_dataset/
│   └── templates/
├── RESULTS/
│   ├── discovery_dataset/
│   ├── pupillometry/
│   ├── replication_dataset/
│   ├── robustness_dataset/
│   └── spin_tests_correlations/
└── FIGURES/                     # Rendered plots, arranged by main/supplementary figure
    ├── main/                    # May be omitted from archives
    └── supplementary/           # May be omitted from archives
```

`CODE/` contains the analysis source. `DATA/` contains processed inputs, metadata, and reference maps. `RESULTS/` contains intermediate Hurst arrays, cleaned tables, statistical outputs, model predictions, and spatial-test results. `FIGURES/` is the destination for exported visualizations; the **figure-generating code** is stored separately under `CODE/figures/`.

## Code: dataset-specific pipelines

### `CODE/datasets/discovery_dataset/`

This pipeline processes the discovery cohort and integrates resting fMRI Hurst estimates with spectroscopy voxel measurements.

**`01_preprocessing/`**

- `1.1_parcellate_fmri.py`: reads externally preprocessed voxelwise fMRI residuals (NIfTI); loads a Schaefer100 atlas; averages signals by parcel; and saves subject-wise MATLAB matrices. Edit the fMRI root, atlas file, session, input filename, and output folder in the script. Preprocessing operations that create the voxelwise residual images occur upstream and are not reproduced here.

**`02_hurst_estimation/`**

- `2.1_estimate_fmri_hurst.m`: discovers subject folders with parcellated `.mat` data, estimates one Hurst value per parcel with the wavelet maximum-likelihood approach, and writes subject-wise MATLAB results. Requires the `nonfractal` and `wmtsa` MATLAB toolboxes plus site-specific input/output paths. **Check the input file selection rule and expected parcellation before running** (see [Known execution caveats](#known-execution-caveats)).

**`03_analyses/`**

- `3.1_mrs_quality_control.ipynb`: cleans MRS measurements, evaluates GABA and Glx spectroscopy-quality criteria (fit error, linewidth/FWHM, signal-to-noise ratio), applies visual-inspection exclusions, and generates `DATA/discovery_dataset/cleaned_MRS_df.csv`.
- `3.2_build_schaefer100_dataframe.ipynb`: consolidates subject-level Hurst `.mat` results into a long-format Schaefer100 table with parcel/network information; creates unfiltered and IQR-filtered CSV outputs.
- `3.3_mrs_voxel_statistics.ipynb`: compares spectroscopy-derived EIB values across MRS voxels using repeated-measures statistical summaries and adjusted pairwise comparisons.
- `3.4_hurst_voxel_statistics.ipynb`: combines the resting Hurst table with the selected cortical parcels belonging to each MRS voxel; computes within-subject voxel averages, descriptives, normality diagnostics, a Friedman omnibus test, and FDR-adjusted paired Wilcoxon comparisons.

**`04_models/`**

- `4.1_von_economo_model.R`: mixed-effects models of resting Hurst, MRS-derived EIB, and von Economo cytoarchitectural classes; exports model predictions and supporting data.
- `4.2_bigbrain_model.R`: models the EIB–Hurst association in relation to a BigBrain-derived cortical gradient, comparing polynomial specifications and exporting predictions.
- `4.3_functional_gradient_model.R`: evaluates linear, quadratic, and cubic relationships between regional Hurst and the principal functional gradient (Margulies Gradient 1).

**`alternative_parcellation_glasser360/`**

- `01_preprocessing/1.1_parcellate_glasser360.py`: analogous fMRI parcellation using Glasser360.
- `02_hurst_estimation/2.1_estimate_glasser360_hurst.m`: estimates regional Hurst values for the alternative parcellation.
- `03_analyses/3.1_build_glasser360_dataframe.ipynb`: merges Glasser360 outputs and creates raw and IQR-cleaned analysis tables.

### `CODE/datasets/replication_dataset/`

The replication workflow analyzes Hurst differences across resting and working-memory task conditions, regional/network structure, and spectroscopy associations. The tasks commonly appear as `Rest`, `0Back`, `1Back`, and `2Back` in downstream tables.

**`01_preprocessing/`**

- `1.1_parcellate_fmri.py`: parcellates preprocessed voxelwise residuals for individual participants and working-memory tasks using Schaefer100. Paths and participant/task identifiers are provided as script configuration and must match local files. As currently written, it loops over 0-, 1-, and 2-back but not resting-state inputs.

**`02_hurst_estimation/`**

- `2.1_estimate_fmri_hurst.m`: MATLAB parcelwise Hurst estimates for task-coded participant files; includes Rest plus three working-memory conditions in its expected input list. It requires the external MATLAB wavelet toolboxes, local parcellated data, and configured output paths.

**`03_analyses/`**

- `3.1_build_schaefer100_dataframe.ipynb`: creates the participant × task × parcel dataset from `.mat` Hurst outputs, attaches parcel/network labels, and applies IQR-based outlier exclusion.
- `3.2_cortical_statistics_schaefer100.ipynb`: cortical task analyses, including parcelwise repeated-measures ANOVAs, two-factor/network-related analyses, network and task contrasts, network-difference summaries, SVD of Hurst patterns, and associated saved tables/arrays.
- `3.3_fsl_one_sample.sh`: FSL-based voxelwise group inference using `fslmerge` and `randomise` with threshold-free cluster enhancement (TFCE). **This is a separate voxelwise pathway**, requires NIfTI images and a mask not included here, and has local paths that must be edited.

**`04_models/`**

- `4.1_neurometabolites_models.R`: mixed-effects associations between task-dependent Hurst, EIB, GABA, and Glx from the spectroscopy/fMRS-derived input table; exports model data, fixed effects, and covariance matrices for downstream plots.
- `4.2_functional_gradient_models.R`: compares linear/quadratic/cubic relationships of parcel Hurst with the principal cortical functional gradient.
- `4.2_functional_gradient_models_with_interaction.R`: extends the gradient models with task interactions; supports evaluation of load-dependent gradient–Hurst associations.

**`alternative_parcellation_glasser360/`**

- `01_preprocessing/1.1_parcellate_glasser360.py`: reproduces the parcellation step using the HCP-MMP/Glasser360 atlas.
- `02_hurst_estimation/2.1_estimate_glasser360_hurst.m`: reruns regional Hurst estimation on Glasser360 time series.
- `03_analyses/3.1_build_glasser360_dataframe.ipynb`: creates the alternative-atlas Hurst tables.
- `03_analyses/3.2_cortical_statistics_glasser360.ipynb`: alternative-atlas parcelwise tests and decomposition analyses (including ANOVA and SVD) for parcellation sensitivity.

### `CODE/datasets/robustness_dataset/`

This component examines an additional resting-state cohort (HCP) from an existing MATLAB matrix of Hurst estimates; **there is no separate raw fMRI preprocessing/estimation stage in this folder**.

- `03_analyses/3.1_build_schaefer100_rest_dataframe.ipynb`: converts `DATA/robustness_dataset/Hurst_rest_all.mat` to a parcelwise long-format table and applies IQR cleaning.
- `04_models/4.1_functional_gradient_model.R`: fits polynomial mixed-effects models of resting Hurst versus the principal functional gradient and exports predictions.

## Code: eye tracking

### `CODE/eye_tracking/01_preprocessing/`

The pupil and microsaccade workflows use EyeLink `.edf` data. These recordings and task event files are **not** included.

- `1.1_preprocess_pupil_size.py`: reads EyeLink recordings; interpolates blink periods; identifies velocity-related artifacts; applies temporal filtering; and stores cleaned pupil-size measurements. Input and output directories are site-specific.
- `1.2_extract_microsaccades.py`: reads saccade events from EyeLink recordings and exports timing and eye-position features used in the microsaccade analyses. The script's event selection follows the original analysis convention; review operational definitions before applying it to new data.
- `1.3_residualize_eye_signals.ipynb`: removes modeled task-evoked responses from pupil-size and microsaccade time series, creates residual signals, and prepares analysis-specific temporal sampling (including the `decimated_2s` products).

### `CODE/eye_tracking/02_hurst_estimation/`

- `2.1_estimate_eye_hurst.m`: estimates Hurst exponents from prepared pupil or microsaccade signals using `bfn_mfin_ml`. Choose the appropriate input directory and output directory for each modality; this is a configurable script rather than an automatically generated batch for both signals.

### `CODE/eye_tracking/03_analyses/`

- `3.1_build_eye_dataframe.ipynb`: transforms saved MATLAB Hurst matrices into long-format pupil/microsaccade CSV dataframes with task, modality, participant, and session-order metadata.

### `CODE/eye_tracking/04_models/`

- `4.1_session_order.R`: mixed-effects analysis of session/order effects in eye-tracking Hurst estimates, with model summaries and data for plotting.
- `4.2_peripheral_ei.R`: peripheral Hurst associations with spectroscopy-derived EIB, GABA, and Glx, including multiple-testing adjustments and saved model outputs.

The relevant saved tables and model files are organized under `RESULTS/pupillometry/pupil_size/` and `RESULTS/pupillometry/microsaccades/`. Their exact availability varies between distribution archives.

## Code: figures and shared utilities

### `CODE/figures/main/`

These notebooks assemble and export the main-figure panels from the data and model outputs. They can generally be run independently once their inputs exist; they are not expected to rerun the full statistical workflow.

| Notebook | Figure content |
|---|---|
| `5.1_figure1.ipynb` | Main Figure 1: mean z-scored cortical Hurst maps, network radar plot, spectroscopy-voxel means, EIB/H boxplots, and correlation matrices/scatterplots |
| `5.3_figure2.ipynb` | Main Figure 2: spectroscopy–Hurst models, metabolite/EIB and eye-tracking-related comparisons |
| `5.4_figure3.ipynb` | Main Figure 3: regional cortical Hurst maps, network profiles, task-dependent statistical maps, SVD |
| `5.5_figure4_cytoarchitecture.ipynb` | Main Figure 4: BigBrain and von Economo cytoarchitecture panels |
| `5.6_figure4_functional_gradient.ipynb` | Main Figure 4: functional-gradient fits across discovery, replication, and robustness cohorts |

The numbering is retained from the analysis project (for example, no `5.2_...` notebook is present in this archive); the notebook filenames are not a universal execution order.

### `CODE/figures/supplementary/`

| Notebook | Supplementary content |
|---|---|
| `5.1_s01_gradient_distribution.ipynb` | S1: cortical gradient distributions and cytoarchitectural classes |
| `5.2_s02_voxel_correlations.ipynb` | S2: spectroscopy–Hurst relationships within MRS voxels |
| `5.3_s04_neurometabolites_and_eye.ipynb` | S4: metabolite/EIB and peripheral-signal analyses |
| `5.4_s05_spatial_tests_schaefer100.ipynb` | S5: Schaefer100 spatial/group statistical tests |
| `5.5_s06_pairwise_schaefer100.ipynb` | S6: pairwise cortical comparisons under Schaefer100 |
| `5.6_s07_network_comparisons_schaefer100.ipynb` | S7: network-level session and task contrasts |
| `5.7_s08_svd_gradient_schaefer100.ipynb` | S8: SVD and principal-gradient relations |
| `5.8_s9_discovery_glasser360.ipynb` | S9: discovery dataset, Glasser360 sensitivity analysis |
| `5.9_s10_replication_glasser360.ipynb` | S10: replication dataset, Glasser360 sensitivity analysis |
| `5.10_s11_mrs_cytoarchitecture.ipynb` | S11: spectroscopy voxel correspondence with cortical microstructure |

The S3 label is not represented by a standalone notebook in the supplied code. For final manuscript submission, check the panel numbering and captions against the definitive manuscript.

### `CODE/utils/`

The `utils` package centralizes functionality used by notebooks and scripts. It should be importable by adding the repository's `CODE/` folder to the Python import path (many notebooks do this automatically).

| Module | Main purpose |
|---|---|
| `__init__.py` | Marks `utils` as an importable Python package |
| `helper_functions.py` | NIfTI-to-parcel parcellation and IQR-based outlier filtering |
| `statistics.py` | Parcelwise summaries, network/task statistics, nonparametric comparisons and multiple-testing helpers |
| `mrs.py` | Spectroscopy-value parsing and within-subject comparison utilities |
| `model_helpers.py` | Extracts task-specific intercepts and slopes for fitted-model visualization |
| `network_profiles.py` | Maps parcels to network labels; aggregates cortical/network profiles and decomposition results |
| `network_plots.py` | Network profile and radar-plot functions |
| `surface.py` | Mesh/annotation utilities and cortical surface plotting |
| `svd.py` | Builds subject/task × parcel matrices, processes missing entries, calculates and stores SVD outputs |
| `contrasts.py` | Parses pairwise condition contrasts and prepares contrast plotting matrices |

**Notebooks should be opened from within the repository:** the directory structure matters because scripts identify the project root using the presence of `CODE/` and `DATA/`.

## Data, results, and figures

### `DATA/`: processed scientific inputs and reference information

| Folder / file | Description |
|---|---|
| `DATA/discovery_dataset/20260529_101516_WAND_ses-05_GANNET_MRS_results.csv` | Discovery spectroscopy results before the repository's MRS quality-control step; includes metabolite values and QC fields |
| `DATA/discovery_dataset/exclusions.csv` | Additional MRS exclusions, including manually reviewed records |
| `DATA/discovery_dataset/cleaned_MRS_df.csv` | QC-processed discovery spectroscopy table, read by downstream statistics/models |
| `DATA/discovery_dataset/selected_parcels_overlap30.csv` | Correspondence of MRS voxels to Schaefer100 parcels, including overlap metadata and a 30% selection flag |
| `DATA/replication_dataset/EIB_12s_MRS_static_&_AUC_New_criteria.csv` | Replication spectroscopy/fMRS-derived GABA, Glx and EIB values and inclusion flags, indexed by subject/run |
| `DATA/replication_dataset/selected_parcels_overlap30.csv` | MRS voxel-to-parcel matching for the replication cohort (DLPFC) |
| `DATA/robustness_dataset/Hurst_rest_all.mat` | Derived resting-state Hurst estimates for the HCP robustness cohort |
| `DATA/cytoarchitecture/parcel_cyto_df.csv` | Parcelwise von Economo cytoarchitecture class mapping |
| `DATA/cytoarchitecture/bb_gradient_df.csv` | Parcelwise BigBrain-derived microstructural gradient (the value column retains the historical name `Cyto_class`) |
| `DATA/gradients/Margulies_gradient1_Schaefer100_7Networks.csv` | First functional-gradient value for Schaefer100 parcels |
| `DATA/gradients/Margulies_gradients_Schaefer100_7Networks.csv` | Multiple functional gradients and corresponding parcel labels |
| `DATA/pupillometry/` | Eye-tracking inputs, event regressors, preprocessed signals and/or residual time series as required by the relevant notebooks; large files may be omitted |
| `DATA/templates/` | External cortical atlas volumes, surface annotations, atlas centroids, parcel labels and Glasser-to-network mapping files; can be empty in a size-reduced archive |

**Important:** `DATA/templates/` is a required dependency for several dataframe-building and surface-visualization notebooks. Atlas resources must correspond to the same coordinate space, parcel indexing, and mesh convention as the input data. Merely labeling two files “MNI” does not guarantee identical voxel grids; verify spatial alignment and use appropriate resampling where needed. Do not resample continuous fMRI signals and discrete atlas labels using the same interpolation method by default.

Examples of template filenames referenced in the code include:

```text
Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv
Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz
lh.Schaefer2018_100Parcels_7Networks_order.annot
rh.Schaefer2018_100Parcels_7Networks_order.annot
Glasser360_2mm.Centroid_RAS.csv
Glasser360_Yeo7.csv
lh.HCP-MMP1.annot
rh.HCP-MMP1.annot
```

These filenames are provided to identify dependencies; this repository does not guarantee that the corresponding resources are included or redistributable. The original parcellation-specific atlas provenance and license should be checked before adding the files to a public release.

### `RESULTS/`: Hurst estimates, statistics and exported model outputs

#### `RESULTS/discovery_dataset/`

- `hurst/raw/parcelwise/FIN_D4_modwt_1_5_padding_0_decimate_0/`: subject-level `.mat` Hurst outputs generated from Schaefer100 time series. The term `raw` here means **unaggregated Hurst estimates**, not raw MRI data.
- `hurst/dataframes/`: full and IQR-cleaned Schaefer100 long-format Hurst tables.
- `hurst_glasser360/raw/` and `hurst_glasser360/dataframes/`: equivalent Glasser360 Hurst estimates and tables.
- `statistics/EIB_repeated_measures/`: MRS EIB long/wide tables, omnibus and FDR-adjusted paired tests.
- `statistics/Hurst_repeated_measures/`: voxel-level Hurst summaries, normality diagnostics, Friedman test and paired post-hoc tests.
- `models/`: fitted-model predictions and raw plotting data for von Economo, BigBrain and functional-gradient analyses.

#### `RESULTS/replication_dataset/`

- `hurst/raw/parcelwise/` and `hurst/dataframes/`: Schaefer100 subject/condition regional Hurst estimates and curated long-format tables.
- `hurst_glasser360/raw/parcelwise/` and `hurst_glasser360/dataframes/`: Glasser360 analogues.
- `anova/` and `anova_2way/`: regional repeated-measures/related ANOVA tables.
- `network_within_task/`: network comparisons within each task, with descriptive results and pairwise tests.
- `task_within_network/`: task comparisons within networks, with descriptive results and pairwise tests.
- `network_delta/`: between-network comparisons of task effects, including serialized Python results and exported CSV tables (`csv_tables/`).
- `svd/`: SVD/decomposition products for Schaefer100 (`all/` and `tasks_only/` configurations), including metadata, imputation records, and `.npz` arrays.
- `models/`: spectroscopy/fMRS mixed-effects exports (modeling dataframe, fixed-effect coefficients, variance–covariance matrices for EIB, GABA and Glx) and gradient-model predictions.
- `one_samp_per_session/`: voxelwise group-map outputs for each condition, including stacked four-dimensional NIfTI files and FSL one-sample/TFCE statistical maps.
- `glasser360/`: parallel alternative-atlas statistical outputs (`anova/`, `anova_2way/`, `network_delta/`, `network_within_task/`, `task_within_network/`, and `svd/`, as available).

#### `RESULTS/robustness_dataset/`

- `hurst/dataframes/`: full and IQR-filtered HCP rest Hurst tables.
- `models/`: predictions for the resting Hurst–functional gradient model.

#### `RESULTS/pupillometry/`

- `pupil_size/hurst/raw/decimated_2s/` and `pupil_size/hurst/dataframes/`: Hurst results in `.mat` and `.csv` formats for pupil-size signals.
- `microsaccades/hurst/raw/decimated_2s/` and `microsaccades/hurst/dataframes/`: equivalent outputs for microsaccades.
- Each modality's `models/order/` and `models/MRS/` folders are intended for session-order and metabolite-association outputs. Some of these folders are empty in the supplied archive; the R modeling scripts produce their contents as configured.

#### `RESULTS/spin_tests_correlations/`

- `spin_permutations.npy`: saved spatial rotation/permutation indices.
- `spin_null_distributions.npz`: spatially constrained null distributions.
- `spin_test_results.csv`: resulting spatial correlation test summaries.

Spatial null models help account for the spatial autocorrelation of cortical maps. These saved products support figures and correlation analyses but do not by themselves document every implementation choice; consult the code using them and report parcellation, hemisphere handling, null construction, and permutation count in the manuscript.

### `FIGURES/`: generated figure products

`FIGURES/main/` and `FIGURES/supplementary/` are the intended output directories for saved figure panels and compositions (for example, `Fig_2`, `Fig_3`, `Fig_4`, `Fig_S5`, and `Fig_S10`). The notebooks under `CODE/figures/` define how to reconstruct the plots from their inputs. **A missing figure file in an archive does not imply that its generating code is absent.** If the distributed archive intentionally omits exported PDF/PNG/SVG files for space reasons, retaining the directory structure (e.g., with `.gitkeep` files) can make the intended layout visible in GitHub.

## Installation and dependencies

### Python + R (Conda)

Run these commands in a terminal from the repository root:

```bash
conda env create -f environment.yml
conda activate hurst-neuroimaging
python -m ipykernel install --user --name hurst-neuroimaging --display-name "Python (hurst-neuroimaging)"
```

The provided environment specifies **Python 3.11** and **R 4.4** along with Jupyter, NumPy, pandas, SciPy, Matplotlib, Seaborn, nibabel, Nilearn, scikit-learn, statsmodels, Pingouin, MNE, PyWavelets, BrainSpace, neuromaps, netneurotools, and related packages. It also includes R packages for linear mixed-effects models and data manipulation (`lme4`, `lmerTest`, `emmeans`, `dplyr`, `tidyr`, `readr`, and others). See `environment.yml` for the authoritative dependency list.

For R scripts, use an R installation with the packages in `environment.yml`, and set the working directory to **this repository's root**. In RStudio, check `getwd()` before executing an R script.

> **Dependency reproducibility:** the Conda environment does not pin all package versions and is not a complete lockfile. For a fully reproducible release, export tested package versions and record OS, MATLAB, and FSL versions.

### MATLAB and external estimation packages

MATLAB is needed for the `02_hurst_estimation/` scripts. The code references:

- **nonfractal MATLAB toolbox**, which provides `bfn_mfin_ml`: <https://github.com/wonsang/nonfractal>
- **WMTSA MATLAB toolbox**, used for wavelet time-series analysis: <https://www.atmos.washington.edu/~wmtsa/>
- MATLAB Statistics and Machine Learning functionality, and Parallel Computing Toolbox functionality where `parpool` / `parfor` is enabled.

The MATLAB scripts currently refer to locally installed toolbox directories (e.g., `wmtsa-matlab-0.2.6` and `nonfractal`). Install those tools separately, review their licensing, and edit the `addpath(...)` calls. The external toolboxes are **not** bundled here.

The optional voxelwise permutation path additionally needs **FSL** (`fslmerge`, `randomise`) available on the command line. Eye-tracking preprocessing uses `eyelinkio` to access `.edf` recordings.

## How to reproduce the analyses

### Option A: Start from the supplied intermediate results (recommended)

This is the intended route for readers who do not have access to raw participant data.

1. Download/clone the repository and ensure the **distribution includes the `DATA/` and `RESULTS/` inputs needed for the analysis of interest**. Some results may be distributed as an accompanying archive or deposit instead of tracked GitHub files.
2. Install the Conda environment. Obtain missing cortical templates (`DATA/templates/`) with the same labels, surfaces, parcel order and resolution expected by the notebooks.
3. Work from the repository root. Open an analysis notebook and choose the `hurst-neuroimaging` Python kernel. Existing notebooks search upward for the root containing `CODE/` and `DATA/`; keep those directories in place.
4. Inspect or regenerate the relevant dataframe using a notebook in the dataset's `03_analyses/` folder. **If the published cleaned dataframe already exists, this rebuild is optional.**
5. Run the relevant R model from the repository root to regenerate predictions, coefficient tables and related outputs, if necessary.
6. Open the corresponding notebook in `CODE/figures/main/` or `CODE/figures/supplementary/` and run its cells to save the panels to `FIGURES/`.

For example, to inspect the existing replication Schaefer100 dataframe:

```python
from pathlib import Path
import pandas as pd

root = Path.cwd()  # start Jupyter in the repository root
path = root / "RESULTS/replication_dataset/hurst/dataframes/df_hurst_parcelwise_FIN_D4_modwt_1_5_padding_0_decimate_0_residuals_w_rest_IQR_cleaned.csv"
df = pd.read_csv(path)
print(df[["Subject", "Task", "Parcel", "Network", "Hursts"]].head())
```

To regenerate a model from a compatible R environment:

```bash
# Execute from repository root:
Rscript CODE/datasets/replication_dataset/04_models/4.2_functional_gradient_models.R
```

The figures are **not guaranteed to run without all of their inputs**. For example, surface maps use external atlas annotations; model panels use prediction CSVs; pupil-related panels use eye-tracking model outputs. Reproduce individual analyses in dependency order instead of executing every notebook indiscriminately.

### Option B: Start from your own raw data (methodological code)

The repository also makes the *first stages* of the workflow available even though the original inputs cannot be published. To adapt them:

1. Obtain the appropriate raw or preprocessed NIfTI/EDF data and all ethics/usage permissions from the original data source or your own study.
2. Prepare voxelwise fMRI time series using a preprocessing pipeline compatible with the residual NIfTI file expected by `1.1_parcellate_fmri.py`. The scripts in this repository perform **parcellation of already preprocessed residual fMRI**, not the full motion correction, registration, denoising and confound-regression pipeline from raw scanner data.
3. Choose the correct Schaefer100/Glasser360 atlas, verify coordinate alignment and label correspondence, and edit `atlasfile`, input roots, subject IDs and output roots in the Python parcellation scripts.
4. Run parcellation, inspect the saved `.mat` arrays for the expected `[time points × parcels]` shape, and check the exact MATLAB variable name used by downstream estimators.
5. Install the wavelet estimation toolboxes, configure the corresponding MATLAB script's `data_folder`, `output_folder`, toolbox paths, task list and parameters, then estimate Hurst values.
6. Copy or link Hurst `.mat` outputs into the exact `RESULTS/<dataset>/hurst.../raw/parcelwise/...` location expected by the dataframe-building notebook, or update those notebook input paths consistently.
7. Run downstream analysis notebooks, R models and figure notebooks as in Option A.

For eye tracking, supply compatible EyeLink `.edf` recordings and design/event files; run the pupil preprocessing and saccade extraction scripts, residualize task-driven signals, configure MATLAB Hurst estimation for the chosen signal type, and build the tables.

**Do not assume** that an arbitrary new dataset has compatible sampling rate, task labels, atlas ordering, number of parcels, run duration, or Hurst estimation scale range. Adjust those details scientifically and report any deviations.

### Execution map (inputs and outputs)

| Step | Inputs | Main outputs | Raw recordings needed? |
|---|---|---|---|
| fMRI parcellation | Preprocessed 4D NIfTI; atlas NIfTI | Parcellated `.mat` per participant/condition | Yes, preprocessed source images |
| fMRI Hurst estimation | Parcellated `.mat`; MATLAB toolboxes | Subject/condition parcelwise Hurst `.mat` | No, only parcellated signals (not distributed) |
| Hurst table assembly | Saved Hurst `.mat`; atlas metadata | `RESULTS/.../hurst/dataframes/*.csv` | No |
| Discovery spectroscopy QC | Derived MRS measurements + exclusions | `DATA/discovery_dataset/cleaned_MRS_df.csv` | No, if MRS table provided |
| Cortical/MRS statistical analysis | Cleaned Hurst and spectroscopy tables; parcel correspondence | `RESULTS/.../statistics/`, ANOVA/network/SVD products | No |
| Mixed-effects models | Curated CSV data + gradients/cytoarchitecture | `RESULTS/.../models/*.csv` | No |
| Eye signal preprocessing | EDF + event data | Processed/residual signals | Yes |
| Eye Hurst/dataframe modeling | Residual signals or saved Hurst `.mat`; spectroscopy table | `RESULTS/pupillometry/...` | No if appropriate derived inputs exist |
| Figure generation | Existing tables, models, atlas templates, spatial-test results | `FIGURES/main/`, `FIGURES/supplementary/` | No |

## Data formats and naming conventions

### Typical dataframe columns

| Field | Interpretation |
|---|---|
| `Subject` | Participant identifier (format differs between cohorts) |
| `Task` | Condition label (`Rest`, `0Back`, `1Back`, `2Back` in the replication fMRI tables) |
| `Parcel` | Cortical parcel numeric ID |
| `Network` | Functional network membership for the parcel |
| `Hursts` | Regional fMRI Hurst exponent, as named in the saved fMRI tables |
| `Hurst` | Hurst exponent in eye-tracking dataframe outputs |
| `Likelihood` | Auxiliary value stored in some replication Hurst tables (see estimator/dataframe code for provenance; do not interpret as a universal calibrated probability) |
| `Voxel` | Location or identifier of an MRS voxel |
| `GABA`, `Glx`, `EIB` | Spectroscopy-derived metabolites and E/I proxy, primarily in replication tables |
| `GlxGABA_AlphaTissCorr` | Tissue-corrected discovery spectroscopy Glx/GABA-related measure used for EIB |
| `MarguliesGradient1` | Principal functional-gradient coordinate per parcel |
| `Cyto_class` | Historical column name used by cytoarchitecture reference tables; check the source file to distinguish discrete von Economo class labels from continuous BigBrain values |
| `order` | Session/acquisition ordering variable in peripheral-signal analyses |

Columns and subject-ID conventions are **not identical** across datasets. Use the parsing and merge rules implemented in the appropriate notebook/R file rather than blindly concatenating cohorts.

### Hurst file and method labels

An output path such as

```text
FIN_D4_modwt_1_5_padding_0_decimate_0/
```

encodes the estimator/method label, D4 wavelet filter, MODWT transform, scales 1–5, padding off, and additional decimation off. `IQR_cleaned` denotes a dataframe following the code's interquartile-range outlier filtering. `w_rest` indicates the dataset/export convention including a resting condition where applicable; do not infer that all source datasets include all four tasks.

## Included and missing files

The analyzed ZIP was **intentionally reduced in size**. Accordingly, directory presence and file counts in a particular copy do not define the full research workflow.

- **Source code:** Python, Jupyter, R and MATLAB scripts for the available pipelines are retained under `CODE/`.
- **Raw participant data:** scanner-level fMRI/NIfTI inputs, EyeLink EDF recordings, original event regressors, and other sensitive/raw acquisitions are **not publicly provided**; the corresponding code is retained for methodological transparency.
- **Derived inputs and outputs:** many CSV, `.mat`, `.xlsx`, `.pkl`, and `.npz` files are represented under `DATA/` and `RESULTS/`. Some data and result files may be omitted from a size-limited version. The reader must supply any missing dependencies to run affected stages.
- **Atlas resources:** `DATA/templates/` can be empty in the shared archive, despite downstream notebooks needing its centroids, NIfTI parcellations and cortical surface labels. These are an **external prerequisite**, not automatically created by the Conda installation.
- **Rendered figures:** some or all files under `FIGURES/` may be omitted to limit size. Figure-generating notebooks remain in `CODE/figures/`.
- **Model-derived figure inputs:** some empty `RESULTS/pupillometry/.../models/` folders may require re-running the relevant R scripts, using the necessary processed data, before all figure panels can be regenerated.

Where institutional, legal, or file-size constraints preclude sharing a particular derived input, describe its provenance and access procedure in the accompanying manuscript or dataset deposit. If a reader cannot access that input, the affected analysis is documented but **not fully executable** from the public release alone.

### A note about GitHub and `.gitignore`

The existing `.gitignore` contains broad entries such as `RESULTS/`, `DATA/discovery_dataset/`, `DATA/replication_dataset/`, and `DATA/robustness_dataset/`. These can prevent **untracked** intermediate files from being added to the Git repository, even though they are present in a local ZIP. Review `git status --ignored` and `git check-ignore -v <path>` before committing. To make downstream analyses reproducible from public files, selectively track safe/authorized derived inputs or make them available via an external research-data archive, and document that location here. Do not remove safeguards around sensitive information without reviewing the data first.

## Reproducibility, limitations, and responsible data sharing

### Known execution caveats

The following issues should be reviewed before claiming end-to-end regeneration from the provided scripts:

1. **Discovery parcellation selection:** `CODE/datasets/discovery_dataset/02_hurst_estimation/2.1_estimate_fmri_hurst.m` searches first for filenames containing `Glasser360`, although it is placed in the Schaefer100 workflow, then falls back to any `.mat` file. Confirm that this selects the intended parcellation and exactly one appropriate file per subject.
2. **MATLAB output paths:** the discovery Hurst script has a placeholder output directory; the replication Schaefer100 script writes to `RESULTS/Hurst_results/...`, whereas its dataframe notebook reads `RESULTS/replication_dataset/hurst/raw/...`. Set consistent paths before running. The eye-tracking estimation script also writes to a path different from the one read by `3.1_build_eye_dataframe.ipynb`.
3. **Resting-state task mismatch:** the replication Python parcellation script currently generates only `0back`, `1back`, `2back`, while the subsequent MATLAB script expects the resting `task_00.mat` file too. Supply its preprocessing pathway or adjust the expected list as appropriate.
4. **Atlas dependency:** template and annotation files referenced by notebooks are not necessarily contained in the archive. Verify atlas file versions, cortical parcellation IDs and surface-registration choices.
5. **Case-sensitive directory naming:** some figure notebooks use `hurst_Glasser360` whereas the saved folder is `hurst_glasser360`; this can fail on Linux. Harmonize capitalization before running those figures.
6. **Eye-tracking configuration:** pupil/microsaccade preprocessing and estimation have signal-specific input/output choices and local folders that need to be reconciled; some derived eye-tracking outputs may be absent from the shared archive.
7. **Site-specific paths embedded in source metadata:** some processed CSVs contain original `mask_file` paths. These strings are not portable and may expose workstation information; remove or replace them with portable relative paths in a public release where possible.
8. **Statistical reproducibility:** manuscript methods should record the actual cohort sizes, exclusions, outlier grouping, missing-data handling, effect-size definitions, corrections for multiple tests, model specifications and optimizer diagnostics. A script passing syntax checks alone is not validation of a reported result.

### Suggested quality checks before public release

- Confirm that all files required by each advertised downstream analysis are either committed or available from a documented, permitted data source.
- Validate that all dataframes have the expected keys (`Subject`, `Task`, `Parcel` / `Voxel`), that merges do not unexpectedly duplicate rows, and that Schaefer100/Glasser360 dimensions are correct.
- Compare representative regenerated statistics, mixed-effects coefficients and saved figure panels with the manuscript's reported values.
- Document stochastic seeds, SVD sign choices, spatial-null methods and software versions where these affect reproducibility.
- Exclude `__pycache__/`, `.Rproj.user/`, Jupyter checkpoints, personal absolute paths, participant-identifying metadata and unnecessary intermediate data from the public release.
- Confirm that redistribution is permitted for reference atlases, external software, and all derived participant-level measurements.
- Consider adding a lightweight synthetic-data test of fMRI parcellation, Hurst file parsing and subject/task/parcel join logic.

The purpose of publishing the raw-stage **code**, despite not distributing raw data, is to make methodological decisions inspectable and reusable. The purpose of sharing curated **intermediate results** is to let readers reproduce the parts of the statistical workflow that do not require restricted acquisitions. These are complementary levels of reproducibility and should not be conflated.

## Citation and license

If this repository supports a manuscript, please cite the **final published article** and, where applicable, the software and atlas resources used for analyses (nonfractal/WMTSA, Schaefer atlas, HCP-MMP/Glasser atlas, BigBrain gradient, functional gradients, FSL, and other tools).

**Manuscript citation:** *To be added upon publication.*

**Code license:** No license file was present in the supplied archive. Add an explicit `LICENSE` before public distribution if you intend to grant reuse rights; do not assume that scientific availability automatically specifies permissions. Raw data and third-party atlas/toolbox licensing may differ from the license chosen for the original code.
