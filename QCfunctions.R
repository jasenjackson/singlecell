# Functions required to perform Kamil Slowikowski's QC algorithm for scRNAseq data
# REQUIRED: library(matrixStats) for rowSds()
# REQUIRED: library(scales) for comma()
# REQUIRED: library(ggplot2)
# REQUIRED: library(dplyr)
# REQUIRED: library(MASS)

# calculate basic statistics for each gene
stats_by_gene <- function(count) {
  
  count.mtx <- as.matrix(count)
  gene_stats <- data.frame(
    n = apply(count.mtx, MARGIN=1, function(row) sum(row > 0)),
    mean = rowMeans(count.mtx),
    sd = rowSds(count.mtx)
  ) %>%
    # Calculate the percent of cells/samples that detected each gene
    mutate(percent = 100 * n / ncol(count))
  return (as.data.frame(gene_stats))
}

# calculate basic statistics for each cell
stats_by_cell <- function(mtx, gene_stats, common_genes=TRUE, gene_count=TRUE, umi_count=TRUE, mt_pct=TRUE, verbose=TRUE) {
  
  cell_stats <- data.frame(cell = colnames(mtx))
  
  if (umi_count == TRUE) {
    if (verbose) {cat("Counting the # of reads found in each cell \n")}
    cell_stats$umi_count = apply(mtx, 2, function(col) sum(col))}
  
  if (gene_count == TRUE) {
    if (verbose) {cat("Counting the # of genes expressed in each cell \n")}
    cell_stats$gene_count = apply(mtx, 2,function(col) sum(col > 0))}
  
  if (mt_pct == TRUE) {
    # add "gene" column to input matrix
    if (!("gene" %in% colnames(mtx))) {
      mtx$gene <- rownames(mtx)}
    
    # count # of reads that map to mitochondrial genes
    if (verbose) {cat("Counting the # of mitochondrial reads in each cell \n")}
    cell_stats$mt_count <- apply(X = dplyr::select(filter(mtx, str_detect(gene, "^MT-")), -gene),
                                 MARGIN = 2, function(col) sum(col))
    
    if (verbose) {cat("Calculating % reads mapping to mitochondrial genes for each cell \n")}
    cell_stats <- mutate(cell_stats, percent.mt = (mt_count / umi_count) * 100) %>%
                    dplyr::select(-mt_count)
    mtx <- dplyr::select(mtx, -gene)}
  
  # quantify % of house_keeping (common_genes) expressed in each cell
  if (!(missing(gene_stats)) & (common_genes==TRUE)) {
    if (verbose) {cat("Identifying housekeeping genes (genes expressed in >95% of cells) \n")}
    common_genes <- which(gene_stats$percent > 95) # find common genes
    
    if (verbose) {cat("Calculating % of house keeping genes expressed in each cell \n")}
    cell_stats$percent.cg = apply(
      X = mtx[common_genes,],
      MARGIN = 2,
      FUN = function(col) 100 * sum(col > 0) / length(common_genes))} 
  
  return(cell_stats)
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

# convenience function for density (n=resolution)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# calculate %common genes from original matrix (umi/cpm/log2cpm)
QCstats <- function(count, verbose=TRUE) {
  gene_stats <- stats_by_gene(count)
  cell_stats <- stats_by_cell(count, gene_stats, verbose=verbose)
  return (cell_stats)
}

# plots two geneXumi count scatter plots: percent.cg and density
QCplot <- function(cellstats, hline, vline) {
  # gene x umi count, colored by percent.cg
  p1 <- ggplot(cellstats, aes(umi_count, gene_count, color=percent.cg)) + 
          geom_point(size=0.4) + 
          scale_color_viridis_c(direction = -1)
  
  # gene x umi count, colored by percent.mt
  p2 <- ggplot(cellstats, aes(umi_count, gene_count, color=percent.mt)) + 
    geom_point(size=0.4) + 
    scale_color_viridis_c()
  
  # percent.cg colored by percent.cg
  p3 <- ggplot(cellstats, aes(reorder(cell, percent.cg), percent.cg, color=percent.mt)) +
          geom_point(size=0.4) + 
          scale_color_viridis_c()
  
  # gene x umi count, colored by cell density
  cellstats$density <- get_density(cellstats$umi_count, cellstats$gene_count, n=200)
  p4 <- ggplot(cellstats, aes(umi_count, gene_count, color=density)) + 
    geom_point(size=0.4) + 
    scale_color_viridis_c()
  
  if (!(missing(hline))) {
    p1 <- p1 + geom_hline(yintercept=hline, color="red", size=0.2)
    p2 <- p2 + geom_hline(yintercept=hline, color="red", size=0.2)
    p4 <- p4 + geom_hline(yintercept=hline, color="red", size=0.2)}
  
  if (!(missing(vline))) {
    p1 <- p1 + geom_vline(xintercept=vline, color="red", size=0.2)
    p2 <- p2 + geom_vline(xintercept=vline, color="red", size=0.2)
    p4 <- p4 + geom_vline(xintercept=vline, color="red", size=0.2)}
  
  gridExtra::grid.arrange(p1,p2,p3,p4, ncol=2)
}

# Prepares a cell_stats dataframe for Seurat meta.data slot
QCseurat <- function(cellstats) {
  row.names(cellstats) <- cellstats$cell
  cellstats <- dplyr::select(cellstats, -cell)
  return(cellstats)
}


