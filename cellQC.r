# scRNA cell (and gene) quality control script
# mendillolab.org
# 7-26-2019
# Generates QC plots from gene x cell matrix

# usage: Rscript cellQC.r -m /data/matrix.txt.gz
#        Rscript cellQC.r --matrix /data/matrix.txt.gz

suppressPackageStartupMessages(require(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
# vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf

option_list = list(
  make_option(c("-m", "--matrix"), action="store", default=NA, type="character",
              help="gene count x cell matrix")
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$m==NA) {
  cat("Count matrix: ")
  cat(opt$m)
  cat("\n")
} else {
  cat("Please provide path for a count matrix (-m/--matrix) \n")
}
