library(qvalue)
library(data.table)
library(stringr)
library(plyr)
library(tidyr)
library(tidyverse)
library(readr)
library(readxl)
library(gprofiler2)
library(VennDiagram)
library(igraph)
library(DANDELIONdev)
library(RColorBrewer)

setwd('/project/xuanyao/zining/')

# Extract Secptre Result
trans_discovery_result <- readRDS('./result/sceptre/primary_cd4_t_cell/Stim8hr_sgRNA_asthma_periph/results_run_discovery_analysis.rds')

pval_matrix <- dcast(
  trans_discovery_result[pass_qc == TRUE],
  response_id ~ grna_target,
  value.var = 'p_value'
)

# convert to matrix with rownames
pval_mat <- as.matrix(pval_matrix[, -1])
rownames(pval_mat) <- pval_matrix$response_id

saveRDS(pval_mat, './result/cota_asthma/asthma_periph_200kb_primary_cd4_relevant_genes_trans_assoc.rds')

cat(paste("P-value matrix:", nrow(pval_mat), "response genes x",
          ncol(pval_mat), "gRNA targets\n"))
cat(paste("Non-NA entries:", sum(!is.na(pval_mat)), "\n"))


# WES data
# wes.1 <- fread('./data/wes/UKB_asthma_burden.GCST90085447_buildGRCh38.tsv.gz')
burden_m3 <- fread('./data/wes/UKB_asthma_burden_M3.1.txt')
# Extract gene name
burden_m3[, gene_name := sub("\\(.*", "", Name)]

# Extract p-value as a named vector
wes.p.asthma <- burden_m3 %>% dplyr::select(gene_name, p_value)

# trans QTL data
p.trans <- readRDS('./result/cota_asthma/asthma_periph_200kb_primary_cd4_relevant_genes_trans_assoc.rds')

# remove cross mappable pairs
# cross_mappable_pairs <- fread('../result/cota_ibd/cd_perturb/Nagid_jurkat_essential_genes/Nadig_essential_genes_gwas_overlap_cross_mappable_pairs.tsv')
# for(i in 1:nrow(cross_mappable_pairs)) {
#   gRNA <- cross_mappable_pairs$gRNA_name[i]
#   gene <- cross_mappable_pairs$gene_name[i]
#   
#   # Check if both row and column exist in the matrix
#   if(gene %in% rownames(p.trans) && gRNA %in% colnames(p.trans)) {
#     p.trans[gene, gRNA] <- NA
#   }
# }

p.trans[is.na(p.trans)] <- 1
p.trans[p.trans >= 1] <- 0.99

PVAL_THRESH <- 0.05/nrow(wes.p.asthma)

# Only analyze SNPs with significant GWAS signal
# candidate_gene <- readRDS('./data/perturb_seq/cd4_perturb_seq/Stim8hr_merged/asthma_gwas_burdent_relevant_genes.rds')

genecode <- get(load('./data/gene_position_hg38_v48.rda'))

# Detect Mediating gene
# data(gene_position_hg19, package = "DANDELIONdev")
step1_gene <- med_gene(
  p_trans = as.matrix(p.trans),
  p_burden = wes.p.asthma,
  burden_gene_col = "gene_name",
  burden_p_col = "p_value",
  gene_info = gene_position_hg38,
  gene_id_col = "gene_name",
  periph_candidates = colnames(p.trans),
  periph_mode = "GENE",
  n_cores = 20
)

saveRDS(step1_gene, './result/cota_asthma/asthma_periph_200kb_primary_cd4_essential_genes_med_gene.rds')

# Calculate pairs
step2_gene <- calc_pair(
  mat_sig = step1_gene$mat.sig,
  mat_p = step1_gene$mat.p, p_burden = wes.p,
  periph_candidates = step1_gene$periph_candidate, periph_mode = 'GENE', burden_threshold = PVAL_THRESH,
  burden_gene_col = 'gene_name', burden_p_col = 'p_value'
)


saveRDS(step2_gene, './result/cota_asthma/asthma_periph_200kb_primary_cd4_essential_genes_result_pairs.rds')


step3_gene <- annotation_gene(
  calc_result = step2_gene,
  gene_info = gene_position_hg38,
  gene_id_col = "gene_name"
)

saveRDS(step3_gene, './result/cota_asthma/asthma_periph_200kb_primary_cd4_essential_genes_result_pairs.rds')

chains_gene <- find_gene_chain(step3_gene)

# Network plot

library(data.table)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)


trait_name <- 'Asthma'

step2_gene <- readRDS("/project/xuanyao/zining/result/cota_asthma/asthma_periph_200kb_primary_cd4_essential_genes_result_pairs.rds")

# ── helper: look up gene info from gene_position_hg38 ────────────────────────
get_gene_info <- function(gene_names, ref) {
  ref_dt <- as.data.table(ref)
  lapply(gene_names, function(g) {
    row <- ref_dt[gene_name == g][1]
    if (nrow(row) == 0) return(list(gene_id = NA, info = NA))
    list(
      gene_id = row$gene_ID,
      info    = paste0(row$Chromosome, ':', row$start, '-', row$end,
                       '|', row$type)
    )
  })
}

# ── build result table ────────────────────────────────────────────────────────
dt <- as.data.table(step2_gene$gene_pair)

periph_info <- get_gene_info(dt$periph_cand, gene_position_hg38)
prox_info   <- get_gene_info(dt$med_gene,    gene_position_hg38)

med_assoc <- data.table(
  peripheral_geneid   = sapply(periph_info, `[[`, 'gene_id'),
  peripheral_genename = dt$periph_cand,
  Source              = 'GWAS',            # adjust if you have source info
  peripheral_info     = sapply(periph_info, `[[`, 'info'),
  proximal_geneid     = sapply(prox_info,   `[[`, 'gene_id'),
  proximal_genename   = dt$med_gene,
  proximal_info       = sapply(prox_info,   `[[`, 'info'),
  dandelion_pval      = dt$DANDELION_pval,
  burden_pval         = dt$burden_pval
)


# Strategy 1 ---------------
significance_threshold <- 0.05/nrow(wes.p.asthma)
window_size <- 1e7

# helper: cluster genes by genomic location
cluster_by_location <- function(data, info_col, prefix) {
  data %>%
    mutate(
      chr       = str_extract(.data[[info_col]], "chr[0-9XY]+"),
      start_pos = as.numeric(str_extract(.data[[info_col]], "(?<=:)[0-9]+"))
    ) %>%
    arrange(chr, start_pos) %>%
    group_by(chr) %>%
    mutate(
      dist_to_prev = start_pos - lag(start_pos, default = -1e9),
      is_new_group = dist_to_prev > window_size,
      locus_group  = cumsum(is_new_group)
    ) %>%
    ungroup() %>%
    mutate(locus_id = paste0(prefix, "_", chr, "_Loc_", locus_group))
}

# cluster peripheral and proximal genes
peri_map <- med_assoc %>%
  dplyr::select(peripheral_genename, peripheral_info) %>%
  distinct() %>%
  cluster_by_location("peripheral_info", "Source")

prox_map <- med_assoc %>%
  dplyr::select(proximal_genename, proximal_info) %>%
  distinct() %>%
  cluster_by_location("proximal_info", "Target")

# map back to original data
df_mapped <- med_assoc %>%
  left_join(peri_map, by = c("peripheral_genename", "peripheral_info")) %>%
  dplyr::rename(source_id = locus_id) %>%
  left_join(prox_map, by = c("proximal_genename", "proximal_info"),
            suffix = c('.source', '.target')) %>%
  dplyr::rename(target_id = locus_id)

# build edges
edges_final <- df_mapped %>%
  dplyr::select(from = source_id, to = target_id) %>%
  distinct()

# source nodes (peripheral, red)
nodes_source <- df_mapped %>%
  group_by(source_id) %>%
  summarise(
    label = {
      genes <- sort(unique(peripheral_genename))
      if (length(genes) > 3) {
        paste0(paste(genes[1:3], collapse = "\n"), "\n...")
      } else {
        paste(genes, collapse = "\n")
      }
    },
    out_degree = n_distinct(target_id),
    .groups = "drop"
  ) %>%
  mutate(type = "Peripheral (Locus)", size_val = out_degree * 45) %>%
  dplyr::rename(id = source_id)

# target nodes (proximal, yellow) — label is just gene names joined by newline
nodes_target <- df_mapped %>%
  group_by(target_id) %>%
  summarise(
    label = {
      sub_df <- cur_data() %>%
        dplyr::select(proximal_genename, burden_pval) %>%
        distinct() %>%
        arrange(proximal_genename)
      paste(sub_df$proximal_genename, collapse = "\n")
    },
    in_degree = n_distinct(source_id),
    .groups = "drop"
  ) %>%
  mutate(
    type     = "Proximal (Locus)",
    size_val = ifelse(in_degree <= 5, 0.2, 2 + (in_degree) * 95)
  ) %>%
  dplyr::rename(id = target_id)

# merge nodes and build graph
all_nodes <- bind_rows(nodes_source, nodes_target)
graph     <- tbl_graph(nodes = all_nodes, edges = edges_final)

# ── compute layout first so we can use coordinates for geom_text ──────────────
graph_layout <- create_layout(
  graph,
  layout          = 'graphopt',
  charge          = 0.05,
  mass            = 45,
  spring.length   = 0,
  spring.constant = 1,
  niter           = 30000
)

# ── per-gene significance table with node coordinates ────────────────────────
node_coords <- graph_layout %>%
  filter(type == "Proximal (Locus)") %>%
  dplyr::select(id, x, y)

gene_sig <- df_mapped %>%
  dplyr::select(target_id, proximal_genename, burden_pval) %>%
  distinct() %>%
  dplyr::mutate(significant = !is.na(burden_pval) & burden_pval < significance_threshold) %>%
  left_join(node_coords, by = c("target_id" = "id")) %>%
  group_by(target_id) %>%
  dplyr::mutate(
    # stack gene names vertically below node
    y_label = y - 0.15 - (row_number() - 1) * 0.12
  ) %>%
  ungroup()

# ── plot ──────────────────────────────────────────────────────────────────────
p <- ggraph(graph_layout) +
  
  # edges
  geom_edge_link(
    arrow       = arrow(length = unit(1.4, 'mm'), angle = 15, type = "closed"),
    start_cap   = circle(2.0, 'mm'),
    end_cap     = circle(3, 'mm'),
    color       = "grey60",
    alpha       = 0.7,
    width       = 0.15
  ) +
  
  # nodes
  geom_node_point(aes(color = type, size = size_val)) +
  
  # peripheral (source) node labels — repelled text
  geom_node_text(
    aes(label = label),
    data             = function(x) filter(x, type == "Peripheral (Locus)"),
    repel            = TRUE,
    size             = 2.2,
    color            = "black",
    bg.color         = "white",
    bg.r             = 0.1,
    segment.size     = 0.2,
    max.overlaps     = Inf,
    min.segment.length = 0
  ) +
  
  # proximal (target) non-significant gene labels — black
  geom_text(
    data      = filter(gene_sig, !significant),
    aes(x = x, y = y_label, label = proximal_genename),
    size      = 2.0,
    color     = "black",
    hjust     = 0.5,
    lineheight = 0.9
  ) +
  
  # proximal (target) significant gene labels — red bold
  geom_text(
    data      = filter(gene_sig, significant),
    aes(x = x, y = y_label, label = proximal_genename),
    size      = 2.0,
    color     = "#E41A1C",
    fontface  = "bold",
    hjust     = 0.5,
    lineheight = 0.9
  ) +
  
  # styling
  scale_color_manual(
    values = c("Peripheral (Locus)" = "#A50F15",
               "Proximal (Locus)"   = "#FDBF6F")
  ) +
  scale_size(range = c(0.5, 9), guide = "none") +
  theme_graph() +
  labs(
    title    = paste(trait_name),
    subtitle = "Gene names in RED indicate Burden Test Significance",
    caption  = "Yellow nodes: proximal genes; Red nodes: peripheral GWAS loci.\nSignificant burden test genes shown in red."
  ) +
  theme(legend.position = "bottom")

# p

ggsave(paste0('./result/cota_asthma/asthma_periph_200kb_primary_cd4_essential_genes_network_plot.pdf'), plot = p, width = 15, height = 12, device = cairo_pdf)

saveRDS(df_mapped, file = paste0('./result/cota_asthma/asthma_primary_cd4_essential_genes_network_map_dat.rds'))


