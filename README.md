# Putative genome contamination has minimal impact on the GTDB taxonomy


This repository contains content for an analysis that evaluates the impact of 
contamination on taxonomy. It is not designed for conducting new analyses. 
Instead, it serves as a documented record of the methods used in the manuscript 
and aids in reproducing the results.

For more information, read the manuscript at ... (pending publication).

## 1. Installation

This repository is designed to be run using Luigi on a Unix-based system.

### 1.1 Conda environment

To create the conda environment run the following command to install the required
dependencies specified in `conda/env.yaml`:

```bash
conda env create -f conda/env.yaml
```

### 1.2 Reference data

This analysis requires the GTDB R207 database, including called genes. These
must be downloaded using the published data from here: [https://doi.org/10.48610/85c83e3](https://doi.org/10.48610/85c83e3)

### 1.3 Configuration

The following variables must be set in `workflow/config.py`:

- `DIR_PROJECT`: This is the absolute path to the project directory, it is expected that the reference data are stored in `DIR_PROJECT/output`.

Other variables do not need to be modified assuming the reference data is setup correctly.


## 2. Running the analysis

After you have set the `PYTHONPATH` to the cloned repository, the full analysis can be 
run using the following command:

```bash
python -m workflow run
```

### 2.1 - Running individual tasks

To run individual tasks, you can modify the `__main__.py` file under the `run`
target to specify the output file you want to generate. If you are not familiar
with [Luigi](https://luigi.readthedocs.io/en/stable/api/luigi.task.html), 
you should familiarise yourself with the output target.

### 2.2 - Notable output targets for manuscript

There are a number of key output targets used in generating figures for the 
manuscript, these are detailed below:

**Figure 1:**
- This is an illustrative example.

**Figure 2:**
- 2a: Data are determined by running GUNC (see methods).
- 2b: `./notebook/manuscript/x_n_contigs_y_frequency_hue_pass.ipynb`
- 2c: `./notebook/manuscript/taxonomic_novelty_v2.ipynb`
- 2d: `./notebook/manuscript/gunc_stats_d.ipynb`
- 2e: `./notebook/manuscript/hist_of_top_failed_sp_gunc_stats_e.ipynb`

**Figure 3:**
- This is an illustrative example.

**Figure 4:**
- Directed/random removal: `workflow/fastani_random/e_analyse_results.py`
- Inter-rep distances: `./notebook/manuscript/figure_ani_outcome_rep_rep_dist.ipynb`

**Figure 5:**
- `./workflow/circos_plot/create_circos_plot.py`

**Supplementary tables and main text data**:
- `./workflow/final/v2_generate_master_tsv.py`

**FastANI analysis**:
- `./workflow/fastani_random`

**BAC120 analysis**:
- `./workflow/v2_fasttree_marker_split`

### 2.3 - Notable methods

The majority of the key methods used in this analysis are abstracted into the
`Genome` class, this can be found in `./workflow/model/genome.py`. Note that
these methods are called in their respective Luigi tasks.

**Contig ranking**:
- `get_gunc_max_css_contig_ranking`

**Marker congruence**:
- `get_marker_congruence`

## 3 - Reading individual results

All genome-specific results are stored in the output of the reference data.
To view a specific result, download the reference data and view the file under
each genome directory.
