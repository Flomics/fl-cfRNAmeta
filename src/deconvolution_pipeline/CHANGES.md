# Vorperian nu-SVR deconvolution: changes from upstream

## Origin

The scripts come from the tutorial folder of the original repository:

- Repository: https://github.com/sevahn/deconvolution
- Folder: `deconvolve_cfrna_tutorial/`
- Upstream commit: `5433fe6fd1ee2201a6899f5b2769c1eb1223e614` (2022-02-21)
- Paper: Vorperian et al., *Cell types of origin of the cell-free transcriptome*, Nature Biotechnology 2022
- Basis matrix: `tsp_v1_basisMatrix.txt` (Tabula Sapiens v1, CPM, 62 cell types), not modified

| File | Status |
|---|---|
| `deconvolve.py` | modified (see below) |
| `merge_2.py` | modified (see below) |

The modified `.py` files are kept exactly as they were run on the CRG cluster to produce the
results in `toden/`. That includes a few leftover debug lines, which are listed below as debug-only.
To see the full diff:

```bash
diff -u <upstream>/deconvolve_cfrna_tutorial/deconvolve.py deconvolve.py
diff -u <upstream>/deconvolve_cfrna_tutorial/merge_2.py merge_2.py
```

## `deconvolve.py`

The deconvolution algorithm is unchanged: same scaling, nu-SVR grid
(nu = 0.05, 0.1, 0.15, 0.25, 0.5, 0.75; C = 0.1, 0.5, 0.75, 1, 10), linear kernel, metrics and outputs.
The only change is how the input mixture file is read.

### 1. `processMixture()` reads only the columns it needs

```diff
-def processMixture(mixturePath):
+def processMixture(mixturePath, techrepsrrlist, gene_id='gene_name'):
 ...
-    thuyDat = pd.read_csv(mixturePath, sep = ",", index_col = 0)
+    thuyDat = pd.read_csv(mixturePath, sep = ",", index_col = 0, usecols=[gene_id] + techrepsrrlist)
```

and the call in `main()`:

```diff
-    thuyCPM = processMixture(mixturePath)
+    thuyCPM = processMixture(mixturePath, techrepsrrlist)
```

**Why:** Each deconvolution run handles one sample, but upstream reads the whole sample CPM table
(all genes x all samples) every time. We run all samples in parallel as a SLURM array (243 tasks,
500 MB per task) on a single `sample_cpm.csv` that holds every sample. Reading the full table in
every task wastes memory and I/O and can go over the memory limit. With `usecols`, each task loads only the gene column and its own sample column(s).

**Effect on results:** none. `main()` already kept only `techrepsrrlist` columns right after
reading (`thuyCPM[techrepsrrlist].mean(axis = 1)`), so the data reaching the model is identical.

**Requirement this adds:** the first column (gene names) of the mixture CSV must be named exactly
`gene_name`, or you must pass a different `gene_id`. Upstream accepted any header for the first
column.

### 2. Debug-only (no effect)

- `#pdb.set_trace()` comment added in `main()`.
- Blank lines removed or added around function definitions and at the end of the file.

## `merge_2.py`

The selection logic is unchanged: for each sample, keep the (nu, C) combination with the lowest
RMSE, set negative coefficients to 0, and normalise them to sum to 1.

### 1. Dataset name, output suffix and sample file

```diff
-datasetName = "NCI_cfRNA"
-fend = "_12142021"
-sampData = "samples.csv"
+datasetName = "dec_cfRNAmeta"
+fend = "_results"
+sampData = "sample_cpm.csv"
```

**Why:** These are the user-editable settings at the top of the script. They now point to our
input file and output names. Outputs are `dec_cfRNAmeta_fractions_results.csv` and
`dec_cfRNAmetasupport_vectors__results.csv`. The missing `_` comes from upstream string
concatenation.

### 2. Read only the header of the sample file

```diff
-samps = pd.read_csv(sampData, sep = ",", index_col = 0)
-sampNames = samps.columns
+samps = pd.read_csv(sampData, sep=",", index_col=0, nrows=0)
+sampNames = samps.columns.tolist()
```

**Why:** The script only needs the sample names, which are the column headers. `nrows=0` reads
just the header, so a full CPM matrix with many samples is not loaded into memory. `.tolist()`
turns the pandas Index into a plain list. Iteration and results are the same.

### 3. Debug-only (no effect)

- `os.system("echo DONE")` (4x) and `os.system("echo BIOREP")` inside the loading loop
- `print(samp)` inside the best-hyperparameter loop

These only add progress messages to the SLURM log.

## Known upstream behaviours (not changed, documented for reference)

We left these unchanged on purpose so results stay comparable to the original method:

1. **The CPM threshold is not applied.** In `main()`, `thuyCPM[thuyCPM >= cpmThresh].dropna()`
   is never assigned back to a variable, so `cpmThresh` (we pass 1) has no effect. Only
   genes that are exactly 0 are dropped.
2. **Fit and RMSE are computed in z-scored space.** The basis matrix (per cell type) and the
   mixture are each scaled to mean 0 and variance 1. nu-SVR is fitted on these scaled values. The
   RMSE and Pearson r used to choose (nu, C) compare the scaled mixture with
   `scaledBasis @ normalisedCoefs`.
3. **Hyperparameter selection uses substring matching.** `merge_2.py` collects rows with
   `if samp in i`. A sample name that is a prefix of another (e.g. `S_18` and `S_182`) would pick
   up the other sample's rows. This has been checked for the current `toden/` sample names: there
   are no collisions.
4. **Ties.** If two (nu, C) pairs have exactly the same minimum RMSE, both rows are kept for that
   sample.
5. **Sample names must not contain `-`.** The scripts split names on `-` to recover the sample name.
