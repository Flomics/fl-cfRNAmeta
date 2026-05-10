#!/usr/bin/env Rscript

# src/join_matrices.R
# Joins multiple gene count matrices on 'gene_id' using dplyr.
# Prevents duplicate column headers by dropping them from subsequent files.
# Usage: ./src/join_matrices.R <file1> <file2> [file3 ...] > merged_output.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: join_matrices.R <file1> <file2> [file3 ...] > output.tsv\n", file = stderr())
  quit(status = 1)
}

# Load the first file
first_file <- args[1]
message("Processing: ", first_file)
df <- read_tsv(first_file, show_col_types = FALSE)

if (!"gene_id" %in% colnames(df)) {
  stop("Column 'gene_id' not found in ", first_file)
}

# Join remaining files
for (f in args[-1]) {
  message("Processing: ", f)
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
  
  next_df <- read_tsv(f, show_col_types = FALSE)
  
  if (!"gene_id" %in% colnames(next_df)) {
    stop("Column 'gene_id' not found in ", f)
  }
  
  # Identify columns already present in the merged dataframe (except the join key)
  common_cols <- intersect(colnames(df), colnames(next_df))
  dup_cols <- setdiff(common_cols, "gene_id")
  
  if (length(dup_cols) > 0) {
    message("  Note: The following columns are already present and will be ignored from this file: ", 
            paste(dup_cols, collapse = ", "))
    next_df <- next_df %>% select(-all_of(dup_cols))
  }
  
  # Join on gene_id
  # inner_join ensures we only keep genes present in all batches.
  df <- inner_join(df, next_df, by = "gene_id")
}

message("Finished joining ", length(args), " files.")
message("Final dimensions: ", nrow(df), " genes x ", ncol(df), " columns.")

# Output to stdout
write_tsv(df, stdout())
