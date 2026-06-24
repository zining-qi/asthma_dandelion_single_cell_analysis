library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)
library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)
library(viridis)
library(ggforce)
library(colorspace)
library(ggnewscale)

# Read data
setwd('/project/xuanyao/zining/analysis/')
cell_metadata <- fread('../data/lung_atlas/lung_atlas_cell_metadata.tsv')
asthma_core_genes_df <- fread('../data/asthma/asthma_core_genes_overlapped.csv')
counts <- readMM("../data/lung_atlas/lung_atlas_core_genes_exp_eqtlgen_dgn.mtx")
rownames(counts) <- readLines("../data/lung_atlas/lung_atlas_asthma_core_genes_barcodes.tsv")
colnames(counts) <- asthma_core_genes_df$ensembl_gene_id

umap <- read.table("../data/lung_atlas/lung_atlas_umap_coordinates.tsv", sep = "\t", header = TRUE, row.names = 1)


for(i in 1:nrow(asthma_core_genes_df)) {
  gene_id <- asthma_core_genes_df$ensembl_gene_id[i]  
  gene_name <- asthma_core_genes_df$hgnc_symbol[i]  
  expr <- counts[, gene_id]
  
  
  df <- data.frame(
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    Expression = expr,
    ann_level_1 = cell_metadata$ann_level_1,
    ann_level_3 = cell_metadata$ann_level_3,
    cell_type = paste0(cell_metadata$ann_level_1, ' - ', cell_metadata$ann_level_3)
  )
  
  df$is_expressing <- ifelse(df$Expression > 2, "Yes", "No")
  
  df_filtered <- df %>% 
    filter(!is.na(ann_level_3),
           !ann_level_3 %in% c("", "None", "Unknown", "unknown")) 
  
  #df_filtered$is_expressing <- ifelse(df_filtered$Expression > 2, "Yes", "No")
  
  
  # Define colors
  n_colors <- length(unique(df_filtered$cell_type))
  base_colors <- scales::hue_pal()(n_colors)
  light_colors <- lighten(base_colors, amount = 0)
  transparent_colors <- alpha(light_colors, alpha = 0.5)
  names(transparent_colors) <- sort(unique(df_filtered$cell_type))
  
  p1 <- ggplot(df_filtered, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      aes(fill = cell_type, color = cell_type),
      shape = 21, size = 0.3, stroke = 0  # stroke = 0 removes border
    ) +
    scale_fill_manual(name = "Cell Type", values = transparent_colors) +
    scale_color_manual(values = transparent_colors, guide = "none") + # hide redundant legend
    guides(
      fill = guide_legend(
        override.aes = list(size = 4, shape = 21, alpha = 0.9)  # Explicitly set larger size for legend keys
      )
    )
  # Reset fill/color scales
  p1 <- p1 + new_scale_fill()
  p1 <- p1 + new_scale_color()
  
  expressing_cell_types <- df_filtered %>%
    filter(is_expressing == "Yes") %>%
    pull(cell_type) %>%
    unique()
  
  # Centroids only for those cell types
  # centroids <- df_filtered %>%
  #   filter(cell_type %in% expressing_cell_types) %>%
  #   group_by(cell_type) %>%
  #   summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
  
  centroids <- df_filtered %>%
    filter(is_expressing == "Yes") %>%
    group_by(cell_type) %>%
    group_modify(~ {
      pts <- as.matrix(.x[, c("UMAP_1", "UMAP_2")])
      if (nrow(pts) < 5) {
        return(data.frame(UMAP_1 = median(pts[,1]), UMAP_2 = median(pts[,2])))
      }
      cl <- dbscan(pts, eps = 1.0, minPts = 5)$cluster  # tune eps
      .x %>% mutate(subcluster = cl) %>%
        filter(subcluster > 0) %>%
        group_by(subcluster) %>%
        summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop") %>%
        select(UMAP_1, UMAP_2)
    }) %>%
    ungroup()
  
  # Overlay expressing cells with expression color (also no border)
  # p1 <- p1 +
  #   geom_point(
  #     data = subset(df_filtered, is_expressing == "Yes"),
  #     aes(fill = Expression, color = Expression),
  #     shape = 21, size = 0.3, alpha = 1, stroke = 0
  #   ) +
  #   scale_fill_viridis_c(
  #     option = "magma", direction = -1, begin = 0, end = 1,
  #     name = "Expression", na.value = NA,
  #     guide = guide_colorbar(override.aes = list(shape = 21, size = 3, stroke = 0))
  #   ) +
  #   scale_color_viridis_c(
  #     option = "magma", direction = -1, begin = 0, end = 1,
  #     guide = "none"
  #   ) +
  #   labs(title = gene_name) +
  #   theme_void() +
  #   theme(legend.position = "right")
  
  
  p1 <- ggplot(df_filtered, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      aes(fill = cell_type, color = cell_type),
      shape = 21, size = 0.3, stroke = 0
    ) +
    scale_fill_manual(name = "Cell Type", values = transparent_colors) +
    scale_color_manual(values = transparent_colors, guide = "none") +
    guides(
      fill = guide_legend(
        override.aes = list(size = 4, shape = 21, alpha = 0.9)
      )
    ) +
    new_scale_fill() +
    new_scale_color() +
    # Expression overlay
    geom_point(
      data = subset(df_filtered, is_expressing == "Yes"),
      aes(fill = Expression, color = Expression),
      shape = 21, size = 0.3, alpha = 1, stroke = 0
    ) +
    scale_fill_viridis_c(
      option = "magma", direction = -1, begin = 0, end = 1,
      name = "Expression", na.value = NA,
      guide = guide_colorbar(override.aes = list(shape = 21, size = 3, stroke = 0))
    ) +
    scale_color_viridis_c(
      option = "magma", direction = -1, begin = 0, end = 1,
      guide = "none"
    ) +
    # Labels on top of everything
    geom_text_repel(
      data = centroids,
      aes(x = UMAP_1, y = UMAP_2, label = cell_type),
      size = 2.5,
      color = "black",
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.size = 0.4,
      segment.color = "black",
      segment.curvature = 0,
      force = 2,
      force_pull = 0.5,
      min.segment.length = 0,
      seed = 42
    ) +
    # geom_text_repel(
    #   data = centroids,
    #   aes(x = UMAP_1, y = UMAP_2, label = cell_type),
    #   size = 2.5,
    #   color = "black",
    #   max.overlaps = Inf,
    #   box.padding = 0.3,
    #   segment.size = 0.3,
    #   segment.color = "grey50",
    #   seed = 42
    # ) +
    labs(title = gene_name) +
    theme_void() +
    theme(legend.position = "right")
  
  ggsave(
    filename = paste0('../result/cota_asthma/results/', gene_name, "_umap.png"),
    plot = p1, device = "png", width = 12, height = 6, dpi = 300, bg = "white"
  )
  
  ggsave(
    filename = paste0('../result/cota_asthma/results/', gene_name, "_umap.pdf"),
    plot = p1, device = "pdf", width = 12, height = 6, bg = "white"
  )
  
  
  # Define colors
  # n_colors <- length(unique(df$cell_type))
  # base_colors <- scales::hue_pal()(n_colors)
  # light_colors <- lighten(base_colors, amount = 0)
  # transparent_colors <- alpha(light_colors, alpha = 0.5)
  # names(transparent_colors) <- sort(unique(df$cell_type))
  # 
  # p2 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  #   geom_point(
  #     aes(fill = cell_type, color = cell_type),
  #     shape = 21, size = 0.3, stroke = 0  # stroke = 0 removes border
  #   ) +
  #   scale_fill_manual(name = "Cell Type", values = transparent_colors) +
  #   scale_color_manual(values = transparent_colors, guide = "none") + # hide redundant legend
  #   guides(
  #     fill = guide_legend(
  #       override.aes = list(size = 4, shape = 21, alpha = 0.9)  # Explicitly set larger size for legend keys
  #     )
  #   )
  # # Reset fill/color scales
  # p2 <- p2 + new_scale_fill()
  # p2 <- p2 + new_scale_color()
  # 
  # # Overlay expressing cells with expression color (also no border)
  # p2 <- p2 +
  #   geom_point(
  #     data = subset(df_filtered, is_expressing == "Yes"),
  #     aes(fill = Expression, color = Expression),
  #     shape = 21, size = 0.3, alpha = 1, stroke = 0
  #   ) +
  #   scale_fill_viridis_c(
  #     option = "magma", direction = -1, begin = 0, end = 1,
  #     name = "Expression", na.value = NA,
  #     guide = guide_colorbar(override.aes = list(shape = 21, size = 3, stroke = 0))
  #   ) +
  #   scale_color_viridis_c(
  #     option = "magma", direction = -1, begin = 0, end = 1,
  #     guide = "none"
  #   ) +
  #   labs(title = gene_name) +
  #   theme_void() +
  #   theme(legend.position = "right")
  
  # ggsave(
  #   filename = paste0('../result/cota_asthma/umap_gene_exp_v2/', gene_name, "_umap_without_filter.png"),
  #   plot = p2, device = "png", width = 12, height = 6, dpi = 300, bg = "white"
  # )
  # 
  # ggsave(
  #   filename = paste0('../result/cota_asthma/umap_gene_exp_v2/', gene_name, "_umap_without_filter.pdf"),
  #   plot = p2, device = "pdf", width = 12, height = 6, bg = "white"
  # )
  
  # Print progress
  message("Processed ", gene_name, " (", i, "/", nrow(asthma_core_genes_df), ")")
}



## All genes together

colnames(asthma_core_genes_df) <- c('gene_symbol', 'ensembl_id')

# 1. Prepare expression data for all genes
expr_matrix <- counts[, asthma_core_genes_df$ensembl_id] %>%
  as.matrix() %>%  # Convert sparse to dense matrix
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  pivot_longer(
    cols = -cell_id,
    names_to = "ensembl_id",
    values_to = "Expression"
  ) %>%
  left_join(asthma_core_genes_df, by = "ensembl_id") %>%
  mutate(gene_label = paste0(gene_symbol, "\n(", ensembl_id, ")"))


# 2. Merge with UMAP coordinates and cell types
df_filtered <- df %>% 
  filter(!is.na(ann_level_3),
         !ann_level_3 %in% c("", "None", "Unknown", "unknown"))

plot_data <- df_filtered %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, UMAP_1, UMAP_2, cell_type) %>%
  right_join(expr_matrix, by = "cell_id") %>%
  filter(!is.na(cell_type) & cell_type != "Unknown") %>%
  mutate(is_expressing = ifelse(Expression > 2, "Yes", "No"))

# 3. Define consistent colors
n_colors <- length(unique(plot_data$cell_type))
base_colors <- scales::hue_pal()(n_colors)
light_colors <- lighten(base_colors, amount = 0)
transparent_colors <- alpha(light_colors, alpha = 0.5)
names(transparent_colors) <- sort(unique(plot_data$cell_type))

p1 <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    aes(fill = cell_type, color = cell_type),
    shape = 21, size = 0.3, stroke = 0  # stroke = 0 removes border
  ) +
  scale_fill_manual(name = "Cell Type", values = transparent_colors) +
  scale_color_manual(values = transparent_colors, guide = "none") + # hide redundant legend
  guides(
    fill = guide_legend(
      override.aes = list(size = 4, shape = 21, alpha = 0.9)  # Explicitly set larger size for legend keys
    )
  )
# Reset fill/color scales
p1 <- p1 + new_scale_fill()
p1 <- p1 + new_scale_color()

# Overlay expressing cells with expression color (also no border)
p1 <- p1 +
  geom_point(
    data = subset(plot_data, is_expressing == "Yes"),
    aes(fill = Expression, color = Expression),
    shape = 21, size = 0.3, alpha = 1, stroke = 0
  ) +
  scale_fill_viridis_c(
    option = "magma", direction = -1, begin = 0, end = 1,
    name = "Expression", na.value = NA,
    guide = guide_colorbar(override.aes = list(shape = 21, size = 3, stroke = 0))
  ) +
  scale_color_viridis_c(
    option = "magma", direction = -1, begin = 0, end = 1,
    guide = "none"
  ) +
  labs(title = '') +
  theme_void() +
  theme(legend.position = "right")

# Save the plot
ggsave(
  filename = paste0('../result/cota_asthma/umap_gene_exp_v2/', "all_gene_umap.png"),
  plot = p1, device = "png", width = 12, height = 6, dpi = 300, bg = "white"
)


ggsave(
  filename = paste0('../result/cota_asthma/umap_gene_exp_v2/', 'all_gene_umap.pdf'),
  plot = p1, device = "pdf", width = 12, height = 6, bg = "white"
)


# # 2. Merge with UMAP coordinates and cell types
# plot_data <- df %>%
#   tibble::rownames_to_column("cell_id") %>%
#   dplyr::select(cell_id, UMAP_1, UMAP_2, cell_type) %>%
#   right_join(expr_matrix, by = "cell_id") %>%
#   filter(!is.na(cell_type) & cell_type != "Unknown") %>%
#   mutate(is_expressing = ifelse(Expression > 2, "Yes", "No"))
# 
# # 3. Define consistent colors
# n_colors <- length(unique(plot_data$cell_type))
# base_colors <- scales::hue_pal()(n_colors)
# light_colors <- lighten(base_colors, amount = 0)
# transparent_colors <- alpha(light_colors, alpha = 0.5)
# names(transparent_colors) <- sort(unique(plot_data$cell_type))
# 
# p2 <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
#   geom_point(
#     aes(fill = cell_type, color = cell_type),
#     shape = 21, size = 0.3, stroke = 0  # stroke = 0 removes border
#   ) +
#   scale_fill_manual(name = "Cell Type", values = transparent_colors) +
#   scale_color_manual(values = transparent_colors, guide = "none") + # hide redundant legend
#   guides(
#     fill = guide_legend(
#       override.aes = list(size = 4, shape = 21, alpha = 0.9)  # Explicitly set larger size for legend keys
#     )
#   )
# # Reset fill/color scales
# p2 <- p2 + new_scale_fill()
# p2 <- p2 + new_scale_color()
# 
# # Overlay expressing cells with expression color (also no border)
# p2 <- p2 +
#   geom_point(
#     data = subset(plot_data, is_expressing == "Yes"),
#     aes(fill = Expression, color = Expression),
#     shape = 21, size = 0.3, alpha = 1, stroke = 0
#   ) +
#   scale_fill_viridis_c(
#     option = "magma", direction = -1, begin = 0, end = 1,
#     name = "Expression", na.value = NA,
#     guide = guide_colorbar(override.aes = list(shape = 21, size = 3, stroke = 0))
#   ) +
#   scale_color_viridis_c(
#     option = "magma", direction = -1, begin = 0, end = 1,
#     guide = "none"
#   ) +
#   labs(title = '') +
#   theme_void() +
#   theme(legend.position = "right")
# 
# # Save the plot
# ggsave(
#   filename = paste0('../result/cota_asthma/umap_gene_exp_v2/', "all_gene_umap_without_filter.png"),
#   plot = p2, device = "png", width = 12, height = 6, dpi = 300, bg = "white"
# )
# 
# 
# ggsave(
#   filename = paste0('../result/cota_asthma/umap_gene_exp_v2/', 'all_gene_umap_without_filter.pdf'),
#   plot = p2, device = "pdf", width = 12, height = 6, bg = "white"
# )


