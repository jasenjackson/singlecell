# Functions required to perform Kamil Slowikowski's QC algorithm for scRNAseq data
# REQUIRED: library(matrixStats) for rowSds()
# REQUIRED: library(scales) for comma()

# calculate basic statistics for each gene
stats_by_gene <- function(log2cpm) {
  log2cpm.mtx <- as.matrix(log2cpm)
  gene_stats <- data.frame(
    n = apply(log2cpm.mtx, MARGIN=1, function(row) sum(row > 0)),
    mean = rowMeans(log2cpm.mtx),
    sd = rowSds(log2cpm.mtx)
  ) %>%
    # Calculate the percent of cells/samples that detected each gene
    mutate(percent = 100 * n / ncol(log2cpm))
  return (as.data.frame(gene_stats))
}

# calculate %common genes detected for each sample
# returns transposed df with added column
pcg_from_gene_stats <- function(log2cpm, gene_stats) {

  # Define a set of common genes that are detected in > 95% of samples
  common_genes <- which(gene_stats$percent > 95)

  # For each sample, calculate the percent of common genes detected
  cell_stats <- data.frame(
    cell = colnames(log2cpm),
    percent_common = apply(
      X = log2cpm[common_genes,],
      MARGIN = 2,
      FUN = function(col) 100 * sum(col > 0) / length(common_genes)
    ),
    n_genes_detected = apply(
      X = log2cpm,
      MARGIN = 2,
      FUN = function(col) sum(col > 0)
    )
  )
  # transpose matrix and paste to cell_stats
  log2cpm <- as.data.frame(t(log2cpm))
  cell_stats <- bind_cols(cell_stats, log2cpm)
  return (cell_stats)
}

# plot % common genes detected
plot_percent_common_genes <- function(gene_stats, cell_stats) {
  percent_common <- cell_stats[,1:2]
  common_genes <- which(gene_stats$percent > 95)
  ggplot(percent_common, aes(reorder(cell, percent_common), percent_common)) +
    geom_point() +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    labs(
      y = "Percent", x = "Cells",
      title = sprintf("Percent of %s common genes detected", comma(length(common_genes)))
    )
}

# calculate %common genes from original matrix (umi/cpm/log2cpm)
percent_common_genes <- function(log2cpm) {
  gene_stats <- stats_by_gene(log2cpm)
  cell_stats <- pcg_from_gene_stats(log2cpm, gene_stats)
  out <- list(gene_stats, cell_stats)
  return (out)
}
