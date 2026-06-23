library(AUCell)
library(dplyr)
library(ggplot2)

### Extract AUCell analysis result ###
cells_AUC <- readRDS('../result/cota_asthma/aucell/cells_AUC_only_core.rds')

cells_assignment <- AUCell::AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

cells_assignment$core_genes$aucThr$thresholds
cells_assignment$core_genes$aucThr$selected

core_gene_assigned <- cells_assignment$core_genes$assignment
length(core_gene_assigned)

par(mfrow=c(1,1))

geneSetName <- 'core_genes'
AUCell::AUCell_plotHist(cells_AUC[geneSetName,])
abline(v=0.08, lwd = 2)
text(
  x = 0.08, 
  y = par("usr")[4] * 0.9,  # 90% of top y-axis limit
  labels = sprintf("AUC > %.3f", 0.08),
  pos = 4,  # Position to the left of the point
  col = "black",
  cex = 0.9
)
mtext("AUC Score Histogram by Core", side = 1, line = 4, adj = 0, cex = 0.8)

newSelectedCells <- colnames(cells_AUC@assays@data@listData[["AUC"]])[cells_AUC@assays@data@listData[["AUC"]][geneSetName, ] > 0.08]
length(newSelectedCells)

saveRDS(newSelectedCells, '../result/cota_asthma/aucell/cell_selected_core_genes_0.08.rds')


### Histogram for AUCell analysis result ###

auc_scores <- AUCell::getAUC(cells_AUC)["core_genes", ]
auc_scores <- auc_scores[auc_scores > 0]

ggplot(data.frame(AUC = auc_scores), aes(x = AUC)) +
  geom_histogram(bins = 60, fill = "skyblue", color = "black") +
  labs(title = 'Histogram of AUC for Core Genes',
       y = "Frequency") +
  theme_minimal() 

ggplot(data.frame(AUC = auc_values), aes(x = AUC)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "white") +
  geom_vline(xintercept = cells_assignment[[geneSetName]]$aucThr$selected, 
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = paste("AUC Histogram for", geneSetName),
    x = "AUC Score",
    y = "Number of Cells",
    caption = "AUC Score Histogram by Core Genes"
  ) +
  theme_minimal() +
  theme(
    plot.caption = element_text(hjust = 0, size = 10, color = "gray30"),
    plot.title = element_text(hjust = 0.5)
  )


### Visulization of Cell type with high DPG expression ###

geneSets <- readRDS('/project/xuanyao/zining/data/asthma/asthma_core_pheripheral_genes_core.rds')

# Combine gene names
geneSets$all_genes <- unique(c(geneSets$core, geneSets$peripheral))

# Combine gene IDs
geneSets$all_id <- unique(c(geneSets$core_id, geneSets$peripheral_id))

# If you also want to combine the overlap IDs (though they seem to be subsets)
geneSets$all_id_overlap <- unique(c(geneSets$core_id_overlap, 
                                    geneSets$peripheral_id_overlap))

saveRDS(geneSets, '/project/xuanyao/zining/data/asthma/asthma_core_pheripheral_genes_core.rds')



cell_metadata <- fread('../data/lung_atlas/lung_atlas_cell_metadata_core.tsv')

core_genes_cell <- readRDS('../result/cota_asthma/aucell/cell_selected_core_genes_0.08.rds')

all_genes_cell <- readRDS('../result/cota_asthma/aucell/cell_selected_global_k1_core_peripheral_genes.rds')

core_genes_cell_df <- cell_metadata %>%
  filter(V1 %in% core_genes_cell) %>%
  select(V1, ann_level_1, ann_level_3)


# Count cells by level 3 only
# Convert to data.frame if it's a list
core_genes_cell_df <- as.data.frame(core_genes_cell_df)
core_genes_cell_df$cell_type <- paste0(core_genes_cell_df$ann_level_1, ' - ', core_genes_cell_df$ann_level_3)

# Then count
count_df <- core_genes_cell_df %>% 
  dplyr::count(cell_type) %>% 
  arrange(desc(n))

# Bar plot
ggplot(count_df, aes(x = reorder(cell_type, -n), y = n)) +
  geom_col(fill = "steelblue", width = 0.8) +
  labs(x = "Cell Type", 
       y = "Number of Cells") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0))  # Remove gap below bars




# Calculate proportions and counts
plot_data <- core_genes_cell_df %>%
  dplyr::count(cell_type) %>%
  mutate(
    prop = n / sum(n),  # Calculate proportion
    label = paste0(n)  # Label with count and %
  ) %>%
  arrange(desc(n))  # Sort by frequency

# Create the plot
ggplot(plot_data, aes(x = reorder(cell_type, -n), y = prop)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(
    aes(label = label), 
    vjust = -0.5, 
    size = 3.5,
    lineheight = 0.9  # Better spacing for 2-line label
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.1))  # Space for labels
  ) +
  labs(
    title = "Cell Type Nominated by Core Genes",
    x = "Cell Type",
    y = "Frenquency of Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank()
  )




# # 1. Calculate total cells per cell type in full metadata
# total_cells <- cell_metadata %>%
#   mutate(cell_type = paste0(ann_level_1, " - ", ann_level_3)) %>%
#   dplyr::count(cell_type, name = "total_cells")
# 
# 
# # 2. Calculate cells in all_genes_cell per cell type
# core_genes_cells <- cell_metadata %>%
#   filter(V1 %in% core_genes_cell) %>%
#   mutate(cell_type = paste0(ann_level_1, " - ", ann_level_3)) %>%
#   dplyr::count(cell_type, name = "positive_cells")
# 
# # 3. Merge and calculate proportions
# plot_data <- total_cells %>%
#   left_join(core_genes_cells, by = "cell_type") %>%
#   mutate(
#     positive_cells = replace_na(positive_cells, 0),
#     proportion = positive_cells / total_cells,
#     label = paste0(positive_cells, "/", total_cells, "\n(", scales::percent(proportion, accuracy = 0.1), ")")
#   ) %>%
#   arrange(desc(proportion))
# 
# # 4. Create the plot
# 
# ggplot(plot_data, aes(x = reorder(cell_type, -proportion), y = proportion)) +
#   geom_col(fill = "steelblue", width = 0.7) +
#   scale_y_continuous(
#     labels = scales::percent_format(),
#     expand = expansion(mult = c(0, 0.1))  # Small buffer at top
#   ) +
#   labs(
#     title = "Proportion of Cells in Selected Cells for Each Cell Type",
#     x = "Cell Type",
#     y = "Proportion of Cells"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     panel.grid.major.x = element_blank()
#   )

