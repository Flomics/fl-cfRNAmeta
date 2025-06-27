#!/usr/bin/env Rscript

# read counts matrix
input_file <- commandArgs(trailingOnly = TRUE)[1]
count_matrix <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE, sep = '\t')

column_names <- colnames(count_matrix)
column_indices <- 3:ncol(count_matrix)

# empty data frame to store the results
result_table <- data.frame(Sample = character(),
                           'number_of_genes_contributing_to_1%_of_reads' = numeric(),
                           'number_of_genes_contributing_to_5%_of_reads' = numeric(),
                           'number_of_genes_contributing_to_10%_of_reads' = numeric(),
                           'number_of_genes_contributing_to_50%_of_reads' = numeric(),
                           'number_of_genes_contributing_to_80%_of_reads' = numeric(),
                           stringsAsFactors = FALSE, check.names = FALSE)

percentages <- c(1, 5, 10, 50, 80)

# iterate over each sample column
for (col_index in column_indices) {
  col_name <- column_names[col_index]
  col_values <- count_matrix[, col_index]
  col_values <- col_values[!is.na(col_values)]
  col_values <- col_values[order(col_values, decreasing=TRUE)]

  # calculate
  genes_count <- sapply(percentages, function(percent) {
    threshold <- (percent / 100) * sum(col_values)
    length(which(cumsum(col_values) <= threshold)) + 1
  })

  result_table[nrow(result_table) + 1, ] <- c(col_name, genes_count)
}

write.table(result_table, file = 'genes_contributing_to_percentage_reads.tsv', sep = '\t', row.names = FALSE, quote = FALSE)