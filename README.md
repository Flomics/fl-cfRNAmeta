# fl-cfRNAmeta 

This repository hosts the code used in the analysis performed as part of the manuscript "**Systematic cross-study assessment of RNA-Seq experimental workflows for plasma cell-free transcriptome profiling**" by Tuñí _et al_.


## Repository Structure

```
fl-cfRNAmeta/
├── README.md
├── nextflow/
├── sra_metadata/
├── src/
└── tables/
```

---

## Contents

### 1. `src/` scripts and notebooks

- **`preprocess_metadata_functions.py`**  
  Main script for preprocessing and harmonizing metadata from multiple cfRNA-seq studies.  
  - Loads and standardizes sample-level metadata from various studies.
  - Merges with dataset-level metadata.
  - Applies dataset-specific cleaning, exclusion, and annotation logic.
  - Outputs harmonized per-sample and per-batch metadata tables.

- **`sra_columns_mapping.py`**  
  Helper functions for renaming and standardizing column names and values across metadata tables.

- **`dataset_mappings.json`**  
  JSON file with mappings for dataset names, colours, orders, or values used across scripts.

- **`boxplots_fig2.R`**  
  R script to generate the following boxplots: percentage of spliced reads, percentage of exonic reads, percentage of fragments mapping to the correct gene orientation, NG80, fragment number, percentage of reads mapping to reference human genome, percentage of reads mapping to ERCC spike-ins, effective fragment length distribution, etc. Also used to create the NG80 vs spliced reads scatterplot, and the percent of human reads vs percent of microbial reads scatterplot.

- **`fig_1a_stacked_barplot.R`**  
  R script for creating donor phenotype stacked bar plot.

- **`fig_1b_heatmap.R`**  
  R script for creating the pre-analytical variables heatmap.

- **`fig3c_biotype_stacked_bp.py`**  
  Python script for the Figure 3C RNA biotype stacked bar plot, one bar per dataset batch. Uses
  `barplot_stacked_plot.py`, `dataset_mappings.json` and `gene_biotype_mappings.json`:  
  `python src/fig3c_biotype_stacked_bp.py tables/gene_tpm_norm.tsv tables/sampleinfo_all-batches.tsv tables/gencode_v39_gene_biotypes.tsv figures/fig3c`

- **`perdataset_dimred.sh`**  
  Runs PCA and tSNE on each dataset separately, from the all-batches expression matrix.
  Splits the matrix with `split_matrix_by_dataset.sh` and then calls `dimred_perdataset.py`
  once per dataset, carrying on if one of them fails:  
  `bash src/perdataset_dimred.sh <gene_tpm_norm.tsv> <sampleinfo.tsv> [out_dir]`  
  Plots land in `<out_dir>/plots/<dataset>/{pca,tsne}/`. Datasets sequenced in more than one
  batch are kept together and told apart inside each plot by the `dataset_batch_label` colour;
  set `SPLIT_COL=dataset_batch` to give each batch its own plots instead.

- **`split_matrix_by_dataset.sh`**  
  Splits an all-batches matrix into one matrix and sampleinfo per dataset. Valid on
  `gene_tpm_norm` because it is normalised per sample, so a subset of its columns is exactly
  the TPM of that dataset alone; this does not hold for TMM.

- **`dimred_perdataset.py`**  
  PCA and tSNE of one dataset, coloured by each metadata variable in turn. Dataset colours,
  markers and ordering come from `dataset_mappings.json`, the rest of the palettes from
  `dimred_plots.py`. Variables absent from the sampleinfo are reported and skipped.

- **`dimred_plots.py`**  
  Colour maps, orderings and the scatter used by `dimred_perdataset.py`. Extracted from the
  internal `bioinfo_utils` so this repository runs on its own.

- **`gene_coverage_profile_fig2_tmpH.ipynb`**  
  Jupyter notebook for plotting gene coverage profiles.

- **`join_count_matrix_and_qc_table.ipynb`**  
  Jupyter notebook for merging sliced count matrices and QC tables into a single matrix or QC table file.

- **`merge_fastqs_array_isolate.sh`**  
  Shell script for merging FASTQ files by array or isolate.

- **`ng80.R`**  
  R script to obtain the NG80 metric reported on the manuscript. Needs the count matrix as input file.
  Writes `genes_contributing_to_percentage_reads.tsv` into the current working directory.

- **`make_ng80_tables.sh`**  
  Generates the two NG80 tables in `tables/` from the unfiltered gene-level matrix, chaining
  `filter_rows.py`, `filter_gene_ids.py` and `ng80.R`:  
  `bash src/make_ng80_tables.sh <gene_tpm_norm.tsv> <GENCODEv39_biotype_gene_info.tsv>`

- **`spliced_reads.sh`**  
  Shell script to obtain the number and the % of spliced reads. Needs the deduplicated BAM file as input file.
---

### 2. `nextflow/`

- **Purpose:**  
  Contains configuration files and parameter sets for running nf-core/rnaseq Nextflow pipeline with each dataset.
- **Files:**
  - `base.config`, `base_params.yml`: Base Nextflow configuration and parameters.
  - `smarter.config`, `smarter_v2_params.yml`, `smarter_v3_params.yml`: Configs for SMARTer protocols.
  - `non_smarter.config`, `non_smarter_reverse_params.yml`, `non_smarter_unstranded_params.yml`: Configs for non-SMARTer protocols.
  - `hg38_gencodev39_params.yml`: Parameters for hg38/Gencode v39 reference.
  - `two-color-illumina.config`: Config for two-color Illumina sequencing.

---

### 3. `sra_metadata/`

- **Purpose:**  
  Stores raw and preprocessed metadata files for each study, as well as supplementary tables and GEO series matrix files.
- **Files:**
  - `<dataset>_metadata.csv` / `<dataset>_metadata_preprocessed.csv`: Raw and processed sample metadata.
  - `<dataset>_supp_table_*.xlsx` / `.tsv`: Supplementary tables with additional sample/batch info.
  - `<dataset>_GSE*_series_matrix.txt`: GEO series matrix files for extracting sample annotations.

---

### 4. `tables/`

Contains output and intermediate tables generated by the preprocessing scripts and downstream analyses:

- **`cfRNA-meta_per_sample_metadata.tsv`**  
  Harmonized per-sample metadata table for all included cfRNA-seq datasets.

- **`cfRNA-meta_per_batch_metadata.tsv`**  
  Harmonized per-batch metadata table summarizing batch-level information.

- **`sampleinfo_all-batches.tsv`**  
  Sample information table including all batches.

- **`genes_contributing_to_percentage_reads.tsv`**  
  NG1/5/10/50/80 per sample, over protein-coding genes only. From `src/make_ng80_tables.sh`.

- **`genes_contributing_to_percentage_reads_no_spike_ins.tsv`**  
  Same, over all genes except SIRV/ERCC spike-ins. From `src/make_ng80_tables.sh`.

- **`gencode_v39_gene_biotypes.tsv`**  
  Ensembl (`gene_biotype`) and custom (`custom_biotype`) biotype of every GENCODE v39 gene and
  spike-in, from biomaRt. Input to `src/fig3c_biotype_stacked_bp.py`.

- **`taxa_simple_df_w_batch.tsv`**  
  Simplified taxa table for downstream taxonomic analyses.

#### Updating metadata information in tables

1. Manually change values in `sra_metadata/dataset_metadata.tsv`
2. Re-run `src/preprocess_metadata_functions.py`
3. Commit the updated `cfRNA-meta_per_{batch|sample}_metadata.tsv` files.

---

## License

See LICENSE file for details.

---

## Contact

For questions or contributions, please open an issue.
