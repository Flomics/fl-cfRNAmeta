# Deconvolution plots

One notebook per cell-type deconvolution figure.

| Notebook | Figure |
|---|---|
| `fig4a_celltype_barplot.ipynb` | Fig 4A: stacked barplot of the mean major-cell-type (MCT) fractions per dataset batch |
| `fig4b_chen_heatmap.ipynb` | Fig 4B: clustered heatmap of the MCT fractions in the Chen dataset |
| `suppfig5_mct_boxplot.ipynb` | Suppl Fig 5: boxplots of the MCT fractions per dataset batch |
| `suppfig6_dataset_heatmaps.ipynb` | Suppl Fig 6A–I: one clustered heatmap per dataset |

## How it works

Each notebook:

1. reads the deconvolution results of every dataset and joins them with the sample metadata
   (`tables/cfRNA-meta_per_sample_metadata.tsv`) and the dataset labels (`src/dataset_mappings.json`);
2. selects the major cell types: those reaching at least 10% of the RNA in at least 10% of the
   samples of one dataset;
3. draws the figure with the plotting helpers in `src/` (`barplot_stacked_plot.py`,
   `heatmap_plot.py`, `boxplot.py`).

Steps 1 and 2 are in `deconv_common.py`, shared by the four notebooks.

## Setup

```bash
conda env create -f src/deconvolution_plots/environment.yml
conda activate cfrnameta_deconvolution_plots
```

## Input

Per dataset, the two files written by the nu-SVR deconvolution script:

- `dec_cfRNAmeta_fractions_results.csv`: one row per cell type (62) plus the rows `r` and `rmse`,
  one column per sample
- `sample_list.txt`: the sample IDs, one per line

Datasets: `Flomics_2 chalasani_merged chen ibarra moufarrej reggiardo roskams tao toden zhu decru`.
A dataset whose files are missing is reported and skipped.

Set the paths in the first cell of each notebook. `{dataset}` is replaced by each dataset name:

```python
FRACTIONS_FILE = "/path/to/{dataset}/dec_cfRNAmeta_fractions_results.csv"
SAMPLE_LIST_FILE = "/path/to/{dataset}/sample_list.txt"
OUTPUT_DIR = "../../output/deconvolution"
```

In `suppfig6_dataset_heatmaps.ipynb`, `USE_JULY_2025_METADATA = True` reads the metadata and the
dataset labels from git commit `8598ef0` instead of the current files; set it to `False` to use
the current ones.

## Run

Open the notebook from `src/deconvolution_plots/`, edit the first cell and run all cells.

## Output

Written to `OUTPUT_DIR`, each figure as pdf, png and svg:

| Notebook | Files |
|---|---|
| `fig4a_celltype_barplot.ipynb` | `fig4a_celltype_barplot.*` |
| `fig4b_chen_heatmap.ipynb` | `fig4b_chen_heatmap.*` (heatmap), `fig4b_chen_heatmap_legend.*` (same heatmap with legends) |
| `suppfig5_mct_boxplot.ipynb` | `suppfig5_mct_boxplot.*` |
| `suppfig6_dataset_heatmaps.ipynb` | `suppfig6A_Flomics_2.*` … `suppfig6I_zhu.*` |

`fig4a_celltype_barplot.ipynb` also prints the MCT mean ± SD and the erythrocyte + platelet mean ± SD;
`suppfig5_mct_boxplot.ipynb` prints the Zhu platelet fraction in healthy vs liver-cancer samples with a
Mann–Whitney U test.
