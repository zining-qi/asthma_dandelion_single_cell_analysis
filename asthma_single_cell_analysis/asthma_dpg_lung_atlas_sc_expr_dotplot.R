library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)

#### Dotplot for DPG expression per cell ####

data_dir <- "."  # adjust

counts <- ReadMtx(
  mtx = file.path(data_dir, "matrix.mtx"),
  cells = file.path(data_dir, "barcodes.tsv"),
  features = file.path(data_dir, "features.tsv")
)

seurat_obj <- CreateSeuratObject(counts = counts)

seurat_obj <- NormalizeData(seurat_obj)

meta <- read.csv("metadata.csv", row.names = 1)

meta <- data.frame(
  ann_level_1 = cell_metadata$ann_level_1,
  ann_level_3 = cell_metadata$ann_level_3,
  cell_type = paste0(cell_metadata$ann_level_1, " - ", cell_metadata$ann_level_3)
)

# CRITICAL: set rownames to barcodes
rownames(meta) <- cell_metadata$barcode   # or colnames(seurat_obj) if already aligned

seurat_obj <- AddMetaData(seurat_obj, metadata = meta)
seurat_obj$cell_type <- factor(seurat_obj$cell_type)

# or manual ordering if needed
# seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = desired_order)

DotPlot(
  seurat_obj,
  features = genes,
  group.by = "celltype"   # replace with your column
) + RotatedAxis()


#### Convert h5ad to mtx ####

# h5ad_path <- "./lung_atlas/2b415299-371a-4e95-8cba-8c6036a72ad5.h5ad"
# 
# sce <- readH5AD(h5ad_path)
# assayNames(sce)
# mat <- assay(sce, "X")   # or "counts" if available
# mat <- as(mat, "dgCMatrix")
# 
# gene_names <- rowData(sce)
# asthma_core_genes_df <- fread('../data/asthma/asthma_core_genes_overlapped.csv')
# genes_use <- intersect(asthma_core_genes_df$ensembl_gene_id, gene_names)
# idx <- which(gene_names %in% genes_use)
# 
# mat_sub <- mat[idx, ]
# gene_names_sub <- gene_names[idx]
# 
# gene_names_sub <- gene_names_sub[order(match(gene_names_sub, asthma_genes))]
# mat_sub <- mat_sub[order(match(gene_names_sub, asthma_genes)), ]
# barcodes <- colnames(mat_sub)
# 
# out_dir <- "/project/xuanyao/zining/data/lung_atlas/asthma_subset/"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # matrix
# writeMM(mat_sub, file = paste0(out_dir, "matrix.mtx"))
# 
# # barcodes
# write.table(
#   barcodes,
#   file = paste0(out_dir, "barcodes.tsv"),
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = FALSE
# )
# 
# # features (3-column format required by Seurat)
# features <- data.frame(
#   gene_id = gene_names_sub,
#   gene_name = gene_names_sub,
#   feature_type = "Gene Expression"
# )
# 
# write.table(
#   features,
#   file = paste0(out_dir, "features.tsv"),
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = FALSE
# )


