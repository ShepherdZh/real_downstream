#!/usr/bin/env Rscript

#' CellChat Cell-Cell Communication Analysis Wrapper
#' Analyzes intercellular communication networks from single-cell data

suppressPackageStartupMessages({
  library(optparse)
  library(CellChat)
  library(Seurat)
  library(anndata)
})

# Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input h5ad file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output RDS file for CellChat object", metavar = "character"),
  make_option(c("--species"), type = "character", default = "human",
              help = "Species: 'human' or 'mouse' [default %default]", metavar = "character"),
  make_option(c("--min_cells"), type = "integer", default = 10,
              help = "Minimum cells per cell group [default %default]", metavar = "integer"),
  make_option(c("--cell_type_col"), type = "character", default = "cell_type",
              help = "Column name for cell type annotations [default %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input h5ad file must be specified with --input", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output RDS file must be specified with --output", call. = FALSE)
}

# Read h5ad file
cat("Loading data from", opt$input, "...\n")
adata <- read_h5ad(opt$input)

# Convert to Seurat object (if needed for easier handling)
# Extract expression matrix and metadata
if (opt$cell_type_col %in% colnames(adata$obs)) {
  meta_data <- adata$obs
  cell_labels <- meta_data[[opt$cell_type_col]]
} else {
  stop(paste0("Cell type column '", opt$cell_type_col, "' not found in h5ad metadata.
              Available columns: ", paste(colnames(adata$obs), collapse = ", ")))
}

# Create CellChat object
cat("Creating CellChat object...\n")
data_input <- as.matrix(adata$X)
if (is.null(rownames(data_input))) {
  rownames(data_input) <- adata$var_names
}
if (is.null(colnames(data_input))) {
  colnames(data_input) <- adata$obs_names
}

cellchat <- createCellChat(object = data_input, meta = meta_data, group.by = opt$cell_type_col)

# Load CellChat database
cat("Loading CellChat database for", opt$species, "...\n")
if (opt$species == "human") {
  CellChatDB <- CellChatDB.human
} else if (opt$species == "mouse") {
  CellChatDB <- CellChatDB.mouse
} else {
  stop("Species must be 'human' or 'mouse'")
}

cellchat@DB <- CellChatDB

# Preprocessing
cat("Preprocessing expression data...\n")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project gene expression to protein-protein interaction network
cellchat <- projectData(cellchat, PPI.human)  # or PPI.mouse

# Infer cell-cell communication network
cat("Inferring cell-cell communication network...\n")
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)

# Compute pathway-level communication
cat("Computing communication at pathway level...\n")
cellchat <- computeCommunProbPathway(cellchat)

# Calculate aggregated network
cat("Aggregating cell-cell communication network...\n")
cellchat <- aggregateNet(cellchat)

# Compute centrality measures
cat("Computing network centrality...\n")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Save CellChat object
cat("Saving CellChat object to", opt$output, "...\n")
saveRDS(cellchat, file = opt$output)

cat("CellChat analysis complete!\n")
cat("Summary:\n")
cat("  Cell groups:", length(levels(cellchat@idents)), "\n")
cat("  Significant interactions:", nrow(cellchat@net$count), "\n")
