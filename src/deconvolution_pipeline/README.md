# cfRNA cell-type deconvolution (Vorperian nu-SVR), Flomics version

These scripts estimate how much of a cell-free RNA (cfRNA) sample comes from each cell type. They
use the nu-SVR method from Vorperian et al. with the Tabula Sapiens v1 basis matrix.

This is a lightly modified copy of the tutorial scripts in
[sevahn/deconvolution](https://github.com/sevahn/deconvolution) (`deconvolve_cfrna_tutorial/`,
commit `5433fe6`). The algorithm is unchanged. The changes only affect how input files are read,
so the scripts can run on a large cohort. See [CHANGES.md](CHANGES.md) for each change and why it
was made.

> **Reference:** Vorperian, S. K., Moufarrej, M. N., Tabula Sapiens Consortium & Quake, S. R.
> *Cell types of origin of the cell-free transcriptome.* Nature Biotechnology 40, 855–861 (2022).

## Contents

| File | Purpose |
|---|---|
| `deconvolve.py` | Deconvolves a single sample with nu-SVR over a (nu, C) grid |
| `merge_2.py` | Collects the per-sample results, picks the best (nu, C) for each sample and writes the final fractions |
| `tsp_v1_basisMatrix.txt` | Tabula Sapiens v1 basis matrix (genes x 62 cell types, CPM), unchanged copy of the upstream file |
| `environment.yml` | Tested conda environment for both scripts |
| `CHANGES.md` | Differences from upstream |

## Input requirements

- **Sample table** `sample_cpm.csv`: comma-separated, genes in rows and samples in columns.
  - The first column must be named `gene_name` and contain gene symbols (e.g. `NRGN`), not Ensembl IDs.
  - Values must be CPM and **not log-transformed**.
  - Sample names must not contain `-`.
- **Basis matrix** `tsp_v1_basisMatrix.txt`: included in this folder. The scripts look for it in
  the directory you run them from, so copy `deconvolve.py`, `merge_2.py` and the basis matrix next
  to `sample_cpm.csv` (or run from this folder). To use a different file, edit `BASISMATRIX_PATH`
  at the top of `deconvolve.py`.

## Environment

Create the environment from `environment.yml`:

```bash
conda env create -f environment.yml
conda activate vorperian_deconv
```

It uses the same library versions as upstream (Python 3.7, numpy 1.18.1, pandas 1.0.1,
scikit-learn 0.24.2, scipy 1.5.2, cvxopt 1.2.7), plus scanpy 1.8.1 because `merge_2.py` imports it.
Both scripts were tested in this environment on the upstream tutorial samples and reproduce
upstream's published fractions.

The upstream
[`cfrna_deconv.yml`](https://github.com/sevahn/deconvolution/blob/master/cfrna_deconv.yml) pins
exact package builds from Anaconda's `defaults` channel (e.g. `python=3.7.4=h265db76_1`), so it
may no longer build on every system.

## Usage

Run both steps from the folder that holds `sample_cpm.csv` and the basis matrix.

**1. Deconvolve each sample.** Each call is independent, so you can run them in parallel, for
example as a SLURM job array:

```bash
python3 -c "import deconvolve as deconv; deconv.main(1, 'sample_cpm.csv', ['SAMPLE_ID'], 'nuSVR', 'SAMPLE_ID', jackknife = False)"
```

Arguments: CPM threshold, sample table, the list of columns to average (technical replicates),
the method (`nuSVR`, `QP` or `NNLS`), and the output name.

**2. Merge.** Once every sample has finished, run:

```bash
python3 merge_2.py
```

## How it works

### `deconvolve.py` (one sample)

1. Reads the sample's column from `sample_cpm.csv`, averages any technical replicates and drops
   genes with zero counts.
2. Keeps only the genes the sample shares with the basis matrix and removes duplicate gene names.
3. Scales each basis-matrix column and the sample to mean 0 and variance 1.
4. Fits a linear nu-SVR, with cell types as features and genes as data points, for each of the 30
   combinations of nu ∈ {0.05, 0.1, 0.15, 0.25, 0.5, 0.75} and C ∈ {0.1, 0.5, 0.75, 1, 10}.
5. For each fit, sets negative coefficients to 0, rescales them to sum to 1, rebuilds the sample
   from them, and records the Pearson r and RMSE against the scaled sample.
6. Writes to `<sample>/`:
   - `<sample>_deconvolutionCoefs.csv`: 30 rows, each with the coefficients, `r` and `rmse`
   - `<sample>_coarsegrain___supportVectors.csv`: the genes used as support vectors in each fit

### `merge_2.py` (all samples)

1. Reads the coefficient and support-vector files of every sample listed in the header of
   `sample_cpm.csv`.
2. For each sample, keeps the (nu, C) combination with the lowest RMSE.
3. Sets negative coefficients to 0 and divides by their sum. The result is the fraction of the
   sample's cell-type-specific RNA that comes from each cell type.
4. Writes:
   - `dec_cfRNAmeta_fractions_results.csv`: cell types x samples, with extra `r` and `rmse` rows
   - `dec_cfRNAmetasupport_vectors__results.csv`: support vectors of the selected fit for each sample

To change the output names, edit `datasetName` and `fend` at the top of `merge_2.py`.

## Caveats

These come from upstream and were left unchanged so results stay comparable. More detail is in
[CHANGES.md](CHANGES.md#known-upstream-behaviours-not-changed-documented-for-reference).

- The CPM threshold argument has no effect: only genes with exactly 0 counts are removed.
- Model selection (RMSE) is done on the scaled values, not on the original CPM values.
- `merge_2.py` matches a sample's rows by substring. Make sure no sample name is a prefix of
  another (e.g. `S_18` and `S_182`).
- The fractions describe the relative contribution of **cell-type-specific RNA**, not the
  proportion of cells.
