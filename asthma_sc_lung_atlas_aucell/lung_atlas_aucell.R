library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)
library(data.table)
library(DelayedArray)
library(HDF5Array)
library(DelayedMatrixStats)
library(AUCell)
library(SummarizedExperiment)
library(foreach)
library(doParallel)

basic_path <- '/project/xuanyao/zining/data/'
sc <- readH5AD(paste0(basic_path, "lung_atlas/db3a3799-3c3f-44db-96a5-81da79a5f3e3.h5ad"),
               use_hdf5 = T, reader = 'R')

asthma_genes <- readRDS('/project/xuanyao/zining/data/asthma/asthma_core_peripheral_genes_overlapped_core_atlas.rds')

chunk_size <- 10000
n_cells <- ncol(sc)
rankings_list <- list()

# for (i in seq(1, n_cells, chunk_size)) {
#   chunk <- assay(sc, "X")[, i:min(i+chunk_size-1, n_cells)]
#   chunk_mat <- as(chunk, "dgCMatrix")
#   rankings_list[[i]] <- AUCell_buildRankings(chunk_mat, plotStats = FALSE)
#   print(paste('chunk', i, 'finished'))
# }

chunk_num <- 1
for (i in seq(1, n_cells, chunk_size)) {
  cat(paste0("Processing chunk ", chunk_num, " \n"))
  
  chunk <- assay(sc, "X")[, i:min(i+chunk_size-1, n_cells)]
  chunk_mat <- as(chunk, "dgCMatrix")
  rankings_list[[chunk_num]] <- AUCell_buildRankings(chunk_mat, plotStats = FALSE)
  # Explicitly remove large objects and call garbage collector
  rm(chunk, chunk_mat)
  gc()
  cat(paste0("Chunk ", chunk_num, " finished.\n"))
  chunk_num <- chunk_num + 1
}

rm(sc)             # Remove the object
gc()     

print('Finish ranking for loop.')
# Check which chunks are NULL
#null_chunks <- which(sapply(rankings_list, is.null))
#print(null_chunks)

valid_rankings <- Filter(Negate(is.null), rankings_list)

all_rankings <- lapply(valid_rankings, function(x) assay(x, "ranking"))

combined_rankings <- do.call(cbind, all_rankings)

#stopifnot(ncol(combined_rankings) == sum(sapply(all_rankings, ncol)))

# Step 1: Create a SummarizedExperiment object
se_obj <- SummarizedExperiment(
  assays = list(ranking = combined_rankings),
  colData = DataFrame(cells = colnames(combined_rankings))
)

# Step 2: Add metadata (e.g. nGenesDetected)
metadata(se_obj)$nGenesDetected <- valid_rankings[[1]]@nGenesDetected

# Step 3: Coerce to aucellResults
final_rankings <- as(se_obj, "aucellResults")
ranking_matrix <- getRanking(final_rankings)
#save(ranking_matrix, file = '/project/xuanyao/zining/result/cota_asthma/ranking_mtx_full.RData')

print('Finish merge all ranking.')
  # For parallel processing

# 1. Define chunk size (columns = cells)
n_cells <- ncol(ranking_matrix)
chunk_size <- 10000  # ~11 chunks (adjust as needed)
chunk_indices <- split(1:n_cells, ceiling(seq_along(1:n_cells)/chunk_size))

# 3. Save chunks in parallel (faster)
doParallel::registerDoParallel(cores = 10)  # Use 4 cores

foreach(i = seq_along(chunk_indices)) %dopar% {
   chunk <- ranking_matrix[, chunk_indices[[i]]]
   saveRDS(chunk,  # Save as matrix
           file = paste0("/project/xuanyao/zining/result/cota_asthma/aucell_v3/core_genes/cell_rankings_", i, ".rds"))
}

output_dir <- "/project/xuanyao/zining/result/cota_asthma/aucell_v3/core_genes/cell_ids/"
if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# # 3. Extract and save cell IDs in parallel (with simplified filenames)
foreach(i = seq_along(chunk_indices)) %dopar% {
   # Get cell IDs for the current chunk
  cell_ids <- colnames(ranking_matrix)[chunk_indices[[i]]]

   # Save as TSV with simplified filename "cell_ids_X.tsv"
   write.table(
     data.frame(CellID = cell_ids),
     file = paste0(output_dir, "cell_ids_", i, ".tsv"),  # Removed "_chunk"
     sep = "\t",
     row.names = FALSE,
     quote = FALSE
   )

   # Return success message
   paste("Saved", length(cell_ids), "cell IDs to cell_ids_", i, ".tsv")  # Updated message
}

print('Finish saving rankings.')

# Gene set
geneSets <- list(core_genes = asthma_genes$core_genes_id)

# Claculate AUC
cells_AUC <- AUCell_calcAUC(geneSets, final_rankings)
# save it as
saveRDS(cells_AUC, file = "/project/xuanyao/zining/result/cota_asthma/aucell_v3/cells_AUC_core_genes.rds")

print('Finish AUC calculation.')

# Load libraries
# library(zellkonverter)
# library(SingleCellExperiment)
# library(Matrix)
# library(data.table)
# library(DelayedArray)
# library(HDF5Array)
# library(DelayedMatrixStats)
# library(AUCell)
# library(SummarizedExperiment)
# library(parallel)

# Set up basic path and read data
# basic_path <- '/project/xuanyao/zining/data/'
# sc <- readH5AD(paste0(basic_path, "lung_atlas/2b415299-371a-4e95-8cba-8c6036a72ad5.h5ad"), 
#                use_hdf5 = TRUE, reader = 'R')
# 
# asthma_genes <- readRDS('/project/xuanyao/zining/data/asthma/asthma_core_pheripheral_genes_full.rds')

#expr_matrix <- assay(sc, "X")  # This is likely a DelayedMatrix backed by HDF5

# Convert to sparse matrix once to avoid doing it repeatedly inside the loop
# (optional depending on your memory)
# expr_matrix <- as(expr_matrix, "dgCMatrix")

# Split cell indices into chunks
# chunk_size <- 5000
# n_cells <- ncol(sc)
# chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / chunk_size))
# chunk_numbers <- seq_along(chunks)
# num_cores <- 6

# Define chunk processing function
# process_chunk <- function(cell_idx, chunk_num) {
#   tryCatch({
#     message(paste0("Processing chunk ", chunk_num, " (cells ", min(cell_idx), ":", max(cell_idx), ")"))
#     
#     chunk <- expr_matrix[, cell_idx]
#     chunk_mat <- as(chunk, "dgCMatrix")  # conversion here is OK if small
#     result <- AUCell_buildRankings(chunk_mat, plotStats = FALSE)
#     
#     message(paste0("Completed chunk ", chunk_num, " (", length(cell_idx), " cells)"))
#     return(result)
#   }, error = function(e) {
#     message(paste("Chunk", chunk_num, "failed:", conditionMessage(e)))
#     NULL
#   })
# }

# Run mclapply in parallel (use your desired number of cores)
# rankings_list <- mclapply(
#   chunks, 
#   process_chunk, 
#   mc.cores = num_cores  # adjust this to the number of cores you want to use
# )

# rankings_list <- mclapply(seq_along(chunks), function(i) {
#   process_chunk(chunks[[i]], chunk_num = i)
# }, mc.cores = num_cores)

# Continue as usual
# valid_rankings <- Filter(Negate(is.null), rankings_list)
# all_rankings <- lapply(valid_rankings, function(x) assay(x, "ranking"))
# combined_rankings <- do.call(cbind, all_rankings)
# stopifnot(ncol(combined_rankings) == sum(sapply(all_rankings, ncol)))
# 
# se_obj <- SummarizedExperiment(
#   assays = list(ranking = combined_rankings),
#   colData = DataFrame(cells = colnames(combined_rankings))
# )
# metadata(se_obj)$nGenesDetected <- valid_rankings[[1]]@nGenesDetected
# final_rankings <- as(se_obj, "aucellResults")
# 
# geneSets <- list(
#   core_genes = asthma_genes$core_id_overlap, 
#   peripheral_genes = asthma_genes$peripheral_id_overlap
# )
# 
# cells_AUC <- AUCell_calcAUC(geneSets, final_rankings)
# saveRDS(cells_AUC, file = "../../result/cota_asthma/cells_AUC_full.rds")




