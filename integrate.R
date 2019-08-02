# integrate.R
# 7-26-2019
# Integrates single cell genomics datasets using Seurat methods
# Currently designed to use seurat's reference based CCA method
# for integration. The script expects all datasets to be merged
# into 1 seurat object

# usage: Rscript integrate.r -s "/data/seurat_object.rds" -g "individual" -r "Ind4"
#        Rscript integrate.r --seuratObject" "../data/BreastCancerAtlas.rds" \
#                         --groupLabel "sample_origin" \
#                         --reference "Ind4,Ind5,Ind6,Ind7" \
#                         --out "../data/BreastCancerAtlas.integrated.rds"
library(Seurat)
library(ggplot2)
options(future.globals.maxSize = 4000 * 1024^2)

# parse arguments
suppressPackageStartupMessages(require(optparse))
option_list = list(
  make_option(c("-s", "--seuratObject"), action="store", default=NA, type="character"),
  make_option(c("-g", "--groupLabel"), action="store", default=NA, type="character"),
  make_option(c("-r", "--reference"), action="store", default=NA, type="character"),
  make_option(c("-o", "--out"), action="store", default=NA, type="character")
)
opt = parse_args(OptionParser(option_list=option_list))

main <- function(seuratObject, groupLabel, reference, out) {
  cat("\nSeurat object: "); cat(seuratObject); cat("\n")
  cat("Group label: "); cat(groupLabel); cat("\n")
  cat("Reference(s): "); for (i in reference){ cat(i); cat(" ") }; cat("\n")
  cat("Output path: "); cat(out); cat("\n\n")

  cat("Reading seurat object\n")
  sobj <- readRDS(seuratObject)

  cat("Splitting seurat object by group label\n")
  sobj.list <- SplitObject(sobj, split.by = groupLabel)

  cat("Datasets:");cat(names(sobj.list)); cat("\n")
  cat("Normalizing each dataset\n")
  for (i in names(sobj.list)) {
    cat("\t"); cat(i); cat("\n")
    sobj.list[[i]] <- SCTransform(sobj.list[[i]], verbose = FALSE)
  }

  cat("Preparing for integration\n")
  sobj.features <- SelectIntegrationFeatures(object.list = sobj.list, nfeatures = 3000)
  sobj.list <- PrepSCTIntegration(object.list = sobj.list, anchor.features = sobj.features)

  cat("Setting reference dataset\n")
  reference <- strsplit(reference, split=",")[[1]]
  reference_dataset <- which(names(sobj.list) == reference)

  cat("Finding integration anchors\n")
  sobj.anchors <- FindIntegrationAnchors(object.list = sobj.list, normalization.method = "SCT",
                                         anchor.features = sobj.features, reference = reference_dataset)

  cat("Integrating dataset\n")
  sobj.integrated <- IntegrateData(anchorset = sobj.anchors, normalization.method = "SCT")

  cat("Running Dimension reduction\n")
  sobj.integrated <- RunPCA(object = sobj.integrated, verbose = TRUE)
  sobj.integrated <- RunUMAP(object = sobj.integrated, dims = 1:30)

  cat("Saving integrated object")
  saveRDS(sobj.integrated, file=out)
}

# argument error handling
if (is.na(opt$seuratObject)) {
  cat("Input error: please provide path for a seurat object \n")
} else if (is.na(opt$groupLabel)) {
  cat("Input error: please provide a group label \n")
} else if (is.na(opt$reference)) {
  cat("Input error: please provide the name of your reference dataset \n")
} else if (is.na(opt$out)) {
  cat("Input error: please provide the output path \n")
} else {
  # split up references
  reference <- strsplit(opt$reference, split=",")[[1]]
  main(opt$seuratObject, opt$groupLabel, reference, opt$out)
}
