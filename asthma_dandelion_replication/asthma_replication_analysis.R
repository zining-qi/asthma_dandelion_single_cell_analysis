library(readxl)
library(tidyverse)
library(DANDELIONdev)
library(ggrepel)


asthma_result <- read_excel("result/cota_asthma/table_s2_to_s11_2025-11-12.xlsx", 
                                         sheet = "table s4", skip = 1)

step1_gene <- readRDS('./result/cota_asthma/asthma_periph_200kb_primary_cd4_essential_genes_med_gene.rds')

overlap_periph <- intersect(step1_gene$periph_candidate, unique(asthma_result$`Exposure gene/mapped exposure gene`))

discovery_df <- asthma_result %>% filter(`Exposure gene/mapped exposure gene` %in% overlap_periph)

discovery_core_genes <- unique(discovery_df %>% pull(DPGs)) # 19 core genes

replicate_df <- as.data.frame(step1_gene$mat.p[which(rownames(step1_gene$mat.p) %in% discovery_core_genes), ], keep.rownames = 'response_gene')


# convert to data.table for easy manipulation

# Core gene replication

min_pval  <- apply(replicate_df, 1, min, na.rm = TRUE)


cat(paste("Total core genes tested:", length(min_pval), "\n"))

pval_thresh <- 0.05/nrow(discovery_df)

sum(min_pval < pval_thresh)

# Pair replication

replicate_df_pair_dpgs <- as.data.frame(step1_gene$mat.p[which(rownames(step1_gene$mat.p) %in% discovery_core_genes), ], keep.rownames = 'DPG')

replicate_df_pair_dpgs <- melt(
  replicate_df_pair_dpgs,
  id.vars     = 'DPG',
  variable.name = 'peripheral_gene',
  value.name    = 'p_value'
)

dpg_min <- apply(replicate_df_pair_dpgs, 1, min, na.rm = TRUE)

# ── convert to long format ────────────────────────────────────────────────────
replicate_df_pair <- as.data.table(step1_gene$mat.p, keep.rownames = 'DPG')

replicate_df_pair <- melt(
  replicate_df_pair,
  id.vars     = 'DPG',
  variable.name = 'peripheral_gene',
  value.name    = 'p_value'
)


# all p-values for QQ plot
dpg_min <- data.table(
  DPG      = names(dpg_min),
  min_pval = as.numeric(dpg_min)
)

# all p-values for QQ
all_pvals <- sort(replicate_df_pair$p_value[!is.na(replicate_df_pair$p_value)])
n         <- length(all_pvals)
expected  <- ppoints(n)

qq_dt <- data.table(expected = expected, observed = all_pvals)

# map each DPG min pvalue to its position in the sorted QQ
# map each DPG min pvalue to its position in the sorted QQ
dpg_min[, qq_rank    := match(min_pval, all_pvals)]
dpg_min[, expected_x := expected[qq_rank]]
dpg_min[, observed_y := min_pval]

# ── plot ──────────────────────────────────────────────────────────────────────
# p <- ggplot(qq_dt, aes(x = expected, y = observed)) +
#   geom_point(color = 'grey70', size = 0.4, alpha = 0.5) +
#   geom_abline(slope = 1, intercept = 0,
#               color = 'black', linewidth = 0.6, linetype = 'dashed') +
#   geom_point(data = dpg_min,
#              aes(x = expected_x, y = observed_y),
#              color = '#E41A1C', size = 1.2) +
#   geom_text_repel(
#     data          = dpg_min,
#     aes(x = expected_x, y = observed_y, label = DPG),
#     color         = '#C00000',
#     size          = 3,
#     fontface      = 'bold',
#     max.overlaps  = Inf,
#     segment.size  = 0.3,
#     segment.color = 'grey50'
#   ) +
#   labs(
#     x        = expression(Expected),
#     y        = expression(Observed),
#     title    = 'QQ plot — DANDELION × CD4 perturb-seq replication',
#     subtitle = paste0('n = ', n, ' pairs | ', nrow(dpg_min), ' DPGs labeled')
#   ) +
#   theme_classic(base_size = 12) +
#   theme(plot.title = element_text(face = 'bold'))
# 
# p <- ggplot(qq_dt, aes(x = expected, y = observed)) +
#   geom_point(color = 'grey70', size = 0.4, alpha = 0.5) +
#   geom_abline(slope = 1, intercept = 0,
#               color = 'black', linewidth = 0.6, linetype = 'dashed') +
#   geom_point(data = dpg_min,
#              aes(x = expected_x, y = observed_y),
#              color = '#E41A1C', size = 1.2) +
#   geom_text_repel(
#     data          = dpg_min,
#     aes(x = expected_x, y = observed_y, label = DPG),
#     color         = '#C00000',
#     size          = 3,
#     fontface      = 'bold',
#     max.overlaps  = Inf,
#     segment.size  = 0.3,
#     segment.color = 'grey50'
#   ) +
#   labs(
#     x        = 'Expected p-value (uniform)',
#     y        = 'Observed p-value (DANDELION)',
#     title    = 'QQ plot — DANDELION × CD4 perturb-seq replication',
#     subtitle = paste0('n = ', n, ' pairs | ', nrow(dpg_min), ' DPGs labeled')
#   ) +
#   theme_classic(base_size = 12) +
#   theme(plot.title = element_text(face = 'bold'))
# 
# print(p)
# ggsave('qqplot_dandelion_cd4_replication_original_scale.png', p, width = 6, height = 6, dpi = 300)


# ── rank all DPGs by min pvalue ───────────────────────────────────────────────
# get min pvalue per DPG across all pairs in replicate_df_pair
all_dpg_min <- replicate_df_pair[
  !is.na(p_value),
  .(min_pval = min(p_value)),
  by = DPG
][order(min_pval)]

all_dpg_min[, rank      := .I]
all_dpg_min[, rank_pct  := rank / .N * 100]
all_dpg_min[, is_core   := DPG %in% asthma_core]

n_total <- nrow(all_dpg_min)
cat(paste("Total DPGs ranked:", n_total, "\n"))

# ── rank of 16 asthma core genes ─────────────────────────────────────────────
core_ranks <- all_dpg_min[is_core == TRUE]
setorder(core_ranks, rank)

cat(paste("\nAsthma core genes found in ranking:", nrow(core_ranks), "/ 16\n"))
print(core_ranks[, .(DPG, min_pval, rank, rank_pct)])

# ── how many core genes in top X% ────────────────────────────────────────────
cat("\nCore genes in top percentiles:\n")
for (pct in c(1, 5, 10, 20, 25)) {
  n_in <- sum(core_ranks$rank_pct <= pct)
  cat(paste0("  Top ", pct, "%  (rank <= ", floor(n_total * pct/100), "): ",
             n_in, " / ", nrow(core_ranks), " core genes"))
  if (n_in > 0)
    cat(paste0(" — ", paste(core_ranks$DPG[core_ranks$rank_pct <= pct],
                            collapse=', ')))
  cat("\n")
}

# ── enrichment test: are core genes enriched in the top ranks? ───────────────
for (pct in c(1, 5, 10, 20)) {
  top_n     <- floor(n_total * pct / 100)
  n_core_in <- sum(core_ranks$rank <= top_n)
  n_core    <- nrow(core_ranks)
  # Fisher's exact test
  mat <- matrix(c(n_core_in,
                  top_n - n_core_in,
                  n_core - n_core_in,
                  n_total - top_n - (n_core - n_core_in)),
                nrow = 2)
  ft  <- fisher.test(mat, alternative = 'greater')
  cat(paste0("  Top ", pct, "% enrichment — OR=",
             round(ft$estimate, 2), " p=",
             format(ft$p.value, scientific=TRUE, digits=3), "\n"))
}

p <- ggplot(all_dpg_min, aes(x = rank, y = -log10(min_pval))) +
  geom_point(aes(color = is_core, size = is_core), alpha = 0.6) +
  geom_vline(xintercept = floor(n_total * 0.01),
             linetype = 'dashed', color = 'orange', linewidth = 0.8) +
  geom_vline(xintercept = floor(n_total * 0.005),
             linetype = 'dashed', color = 'red', linewidth = 0.8) +
  annotate("text", x = floor(n_total * 0.01), y = max(-log10(all_dpg_min$min_pval)) * 0.95,
           label = "Top 1%", color = 'orange', hjust = -0.1, size = 3) +
  annotate("text", x = floor(n_total * 0.005), y = max(-log10(all_dpg_min$min_pval)) * 0.85,
           label = "Top 0.5%", color = 'red', hjust = -0.1, size = 3) +
  geom_text_repel(
    data = all_dpg_min[is_core == TRUE],
    aes(label = DPG), size = 2.5, color = '#C00000',
    max.overlaps = Inf, segment.size = 0.2
  ) +
  scale_color_manual(values = c('FALSE' = 'grey70', 'TRUE' = '#E41A1C'),
                     labels = c('Others', 'Asthma DPGs')) +
  scale_size_manual(values  = c('FALSE' = 0.8, 'TRUE' = 2.5), guide = 'none') +
  labs(x = 'Rank (all mediating genes)', y = expression(-log[10](min~p)),
       title = 'Asthma core genes enriched at top ranks',
       subtitle = '15/16 in top 0.5% | 16/16 in top 1%',
       color = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = 'bold'),
        legend.position = 'bottom')

ggsave('dandelion_cd4_replication_rank_plot.png', p, width = 12, height = 8, dpi = 300,device = "png", bg = 'white')

ggsave('dandelion_cd4_replication_rank_plot.pdf', p, width = 12, height = 8, device = "pdf", bg = 'white')

# pcts <- c(0.5, 1, 2, 5, 10, 20, 50)
# enrich_dt <- data.table(
#   percentile = pcts,
#   n_core     = sapply(pcts, function(p) sum(core_ranks$rank_pct <= p)),
#   n_total_in = sapply(pcts, function(p) floor(n_total * p / 100)),
#   expected   = sapply(pcts, function(p) nrow(core_ranks) * p / 100)
# )
# enrich_dt[, observed_pct := n_core / nrow(core_ranks) * 100]
# enrich_dt[, expected_pct := percentile]
# 
# plot_dt <- melt(enrich_dt[, .(percentile, observed_pct, expected_pct)],
#                 id.vars = 'percentile',
#                 variable.name = 'type', value.name = 'pct_core_genes')
# 
# ggplot(plot_dt, aes(x = factor(percentile), y = pct_core_genes, fill = type)) +
#   geom_bar(stat = 'identity', position = 'dodge', width = 0.6) +
#   geom_text(aes(label = paste0(round(pct_core_genes, 1), '%')),
#             position = position_dodge(0.6), vjust = -0.3, size = 3) +
#   scale_fill_manual(values = c('observed_pct' = '#E41A1C',
#                                'expected_pct'  = 'grey70'),
#                     labels = c('Observed', 'Expected (random)')) +
#   labs(x     = 'Top percentile threshold',
#        y     = '% of core genes included',
#        title = 'Enrichment of asthma core genes at top DPG ranks',
#        subtitle = '15/16 in top 0.5% | 16/16 in top 1%',
#        fill  = NULL) +
#   theme_classic(base_size = 12) +
#   theme(plot.title = element_text(face = 'bold'),
#         legend.position = 'bottom')
# 
# set.seed(42)
# 
# # ── sample random genes (same size as core genes, excluding core genes) ────────
# non_core_dpgs <- all_dpg_min$DPG[all_dpg_min$is_core == FALSE]
# n_core        <- nrow(core_ranks)
# 
# # repeat sampling 1000 times to get stable expected distribution
# n_perm  <- 1000
# random_rank_pcts <- replicate(n_perm, {
#   rand_genes <- sample(non_core_dpgs, n_core, replace = FALSE)
#   all_dpg_min[DPG %in% rand_genes, rank_pct]
# })
# 
# # mean rank percentile per random gene position
# random_mean_pct <- rowMeans(random_rank_pcts)
# 
# # ── build plot data ───────────────────────────────────────────────────────────
# # core genes — actual rank percentiles
# core_plot <- data.table(
#   rank_pct = core_ranks$rank_pct,
#   group    = 'Asthma core genes'
# )
# 
# # one random sample for visualization
# rand_sample <- all_dpg_min[DPG %in% sample(non_core_dpgs, n_core)]
# rand_plot   <- data.table(
#   rank_pct = rand_sample$rank_pct,
#   group    = 'Random genes'
# )
# 
# plot_dt <- rbind(core_plot, rand_plot)
# 
# # ── plot 1: boxplot / jitter comparison ──────────────────────────────────────
# p1 <- ggplot(plot_dt, aes(x = group, y = rank_pct, color = group)) +
#   geom_boxplot(width = 0.4, outlier.shape = NA, linewidth = 0.8) +
#   geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
#   geom_text_repel(
#     data = plot_dt[group == 'Asthma core genes'],
#     aes(label = core_ranks$DPG[match(rank_pct, core_ranks$rank_pct)]),
#     size = 2.5, max.overlaps = Inf, segment.size = 0.2
#   ) +
#   scale_color_manual(values = c('Asthma core genes' = '#E41A1C',
#                                 'Random genes'       = 'grey50'),
#                      guide  = 'none') +
#   scale_y_reverse(labels = function(x) paste0(x, '%')) +
#   labs(x     = NULL,
#        y     = 'Rank percentile (lower = better)',
#        title = 'Asthma core genes vs random genes',
#        subtitle = paste0('Core: median rank = ',
#                          round(median(core_ranks$rank_pct), 1),
#                          '% | Random: median = ',
#                          round(median(rand_sample$rank_pct), 1), '%')) +
#   theme_classic(base_size = 12) +
#   theme(plot.title = element_text(face = 'bold'))
# 
# print(p1)
# 
# # ── plot 2: empirical distribution of random median ranks (permutation test) ──
# # compute median rank_pct for each permutation
# random_medians <- apply(random_rank_pcts, 2, median)
# core_median    <- median(core_ranks$rank_pct)
# 
# p_perm <- mean(random_medians <= core_median)
# cat(paste("Permutation p-value (core median rank <=",
#           round(core_median, 2), "%):", p_perm, "\n"))
# 
# p2 <- ggplot(data.table(median_pct = random_medians),
#              aes(x = median_pct)) +
#   geom_histogram(bins = 50, fill = 'grey70', color = 'white') +
#   geom_vline(xintercept = core_median, color = '#E41A1C',
#              linewidth = 1.2, linetype = 'solid') +
#   annotate("text", x = core_median, y = Inf,
#            label = paste0("Core genes\nmedian = ",
#                           round(core_median, 2), "%"),
#            color = '#E41A1C', hjust = -0.1, vjust = 1.5, size = 3.5) +
#   labs(x     = 'Median rank percentile of random gene sets',
#        y     = 'Count',
#        title = 'Permutation test — core gene rank enrichment',
#        subtitle = paste0('Permutation p = ', p_perm,
#                          ' (n = ', n_perm, ' permutations)')) +
#   theme_classic(base_size = 12) +
#   theme(plot.title = element_text(face = 'bold'))
# 
# print(p2)
# 
# ggsave('rank_comparison_core_vs_random.png', p1, width=7, height=6, dpi=150)
# ggsave('rank_permutation_test.png',          p2, width=7, height=5, dpi=150)
# 
# 
# # build summary for bar plot
# bar_dt <- data.table(
#   percentile_label = factor(paste0('Top ', pcts, '%'),
#                             levels = paste0('Top ', pcts, '%')),
#   core_pct         = sapply(pcts, function(p) sum(core_ranks$rank_pct <= p)) / n_core * 100,
#   random_mean_pct  = sapply(seq_along(pcts), function(i) {
#     mean(apply(random_rank_pcts, 2, function(x) sum(x <= pcts[i]))) / n_core * 100
#   }),
#   random_se_pct    = sapply(seq_along(pcts), function(i) {
#     sd(apply(random_rank_pcts, 2, function(x) sum(x <= pcts[i]))) / n_core * 100
#   }),
#   perm_pval        = perm_pvals,
#   pval_label       = ifelse(perm_pvals < 0.001, '***',
#                             ifelse(perm_pvals < 0.01,  '**',
#                                    ifelse(perm_pvals < 0.05,  '*', 'ns')))
# )
# 
# # melt for grouped bar plot
# bar_long <- melt(bar_dt[, .(percentile_label, core_pct, random_mean_pct)],
#                  id.vars      = 'percentile_label',
#                  variable.name = 'group',
#                  value.name    = 'proportion')
# bar_long[, group := ifelse(group == 'core_pct', 'Asthma core genes', 'Random genes')]
# 
# # add error bars for random only
# bar_long[, se := ifelse(group == 'Random genes',
#                         bar_dt$random_se_pct[match(percentile_label,
#                                                    bar_dt$percentile_label)],
#                         NA)]
# 
# p <- ggplot(bar_long, aes(x = percentile_label, y = proportion, fill = group)) +
#   geom_bar(stat = 'identity', position = position_dodge(0.7), width = 0.6) +
#   # error bars for random (mean ± SE)
#   geom_errorbar(
#     aes(ymin = proportion - se, ymax = proportion + se),
#     position = position_dodge(0.7), width = 0.2,
#     na.rm = TRUE, linewidth = 0.6
#   ) +
#   # significance stars above each pair of bars
#   geom_text(
#     data = bar_dt,
#     aes(x = percentile_label, y = pmax(core_pct, random_mean_pct) + 5,
#         label = pval_label),
#     inherit.aes = FALSE,
#     size = 5, color = 'black'
#   ) +
#   scale_fill_manual(values = c('Asthma core genes' = '#E41A1C',
#                                'Random genes'       = 'grey70'),
#                     name = NULL) +
#   scale_y_continuous(limits = c(0, 115),
#                      breaks = seq(0, 100, 25),
#                      labels = function(x) paste0(x, '%')) +
#   labs(x        = 'Top Percentile',
#        y        = '% of genes in top X% rank',
#        title    = 'Asthma DPGs enrichment') +
#   theme_classic(base_size = 12) +
#   theme(plot.title      = element_text(face = 'bold'),
#         legend.position = 'bottom')
# 
# print(p)
# ggsave('rank_enrichment_barplot.png', p, width = 10, height = 6, dpi = 150)
# 
# 
# 
# library(data.table)
# library(ggplot2)
# library(ggrepel)
# 
# discovery_dt <- as.data.table(discovery_df)
# setnames(discovery_dt, 'Exposure gene/mapped exposure gene', 'peripheral_gene')
# 
# # unique data sources
# data_sources <- unique(discovery_dt$data)
# cat(paste("Data sources:", paste(data_sources, collapse=', '), "\n"))
# 
# # ── for each data source, rank core genes by min trans_p ──────────────────────
# results_list <- list()
# 
# for (src in data_sources) {
#   src_dt <- discovery_dt[data == src]
#   
#   # rank core genes by min trans_p
#   core_rank_src <- src_dt[DPGs %in% discovery_core_genes,
#                           .(min_pval = min(trans_p, na.rm=TRUE)),
#                           by = DPGs][order(min_pval)]
#   core_rank_src[, rank_in_src   := .I]
#   core_rank_src[, n_total_src   := nrow(core_rank_src)]
#   core_rank_src[, rank_pct_src  := rank_in_src / n_total_src * 100]
#   core_rank_src[, data_source   := src]
#   
#   results_list[[src]] <- core_rank_src
# }
# 
# src_ranks <- rbindlist(results_list)
# cat("\nCore gene ranks per data source:\n")
# print(src_ranks)
# 
# # ── perturb-seq ranking (from previous analysis) ──────────────────────────────
# # all_dpg_min has rank_pct from CD4 perturb-seq
# perturb_ranks <- all_dpg_min[is_core == TRUE,
#                              .(DPGs = DPG, min_pval_perturb = min_pval,
#                                rank_pct_perturb = rank_pct)]
# 
# # ── merge discovery ranks with perturb-seq ranks ──────────────────────────────
# combined <- merge(src_ranks, perturb_ranks, by = 'DPGs', all.x = TRUE)
# cat("\nCombined ranks:\n")
# print(combined[, .(DPGs, data_source, rank_pct_src, rank_pct_perturb, min_pval, min_pval_perturb)])
# 
# # ── random genes comparison ───────────────────────────────────────────────────
# set.seed(42)
# n_core  <- length(discovery_core_genes)
# n_perm  <- 1000
# 
# # for each data source, compute rank correlation for random gene sets
# rank_cor_list <- list()
# 
# # use all non-core genes available in perturb-seq as the random pool
# non_core_perturb <- all_dpg_min$DPG[!all_dpg_min$DPG %in% discovery_core_genes]
# 
# for (src in data_sources) {
#   src_dt  <- discovery_dt[data == src]
#   core_src <- combined[data_source == src & !is.na(rank_pct_perturb)]
#   n_core_src <- nrow(core_src)
#   
#   obs_cor <- cor(core_src$rank_pct_src, core_src$rank_pct_perturb,
#                  method = 'spearman', use = 'complete.obs')
#   
#   rand_cors <- replicate(n_perm, {
#     # sample random genes from perturb-seq pool
#     rand_genes   <- sample(non_core_perturb, n_core_src, replace = FALSE)
#     rand_perturb <- all_dpg_min[DPG %in% rand_genes,
#                                 .(DPGs = DPG,
#                                   rank_pct_perturb = rank_pct)][order(rank_pct_perturb)]
#     rand_perturb[, rank_pct_src := seq(100/n_core_src, 100, length.out=.N)]
#     if (nrow(rand_perturb) < 3) return(NA)
#     cor(rand_perturb$rank_pct_src, rand_perturb$rank_pct_perturb,
#         method = 'spearman', use = 'complete.obs')
#   })
#   
#   perm_p <- mean(rand_cors <= obs_cor, na.rm = TRUE)
#   cat(paste0(src, ": r=", round(obs_cor,3),
#              " | rand mean=", round(mean(rand_cors,na.rm=TRUE),3),
#              " | p=", round(perm_p,3), "\n"))
#   
#   rank_cor_list[[src]] <- data.table(
#     data_source = src,
#     obs_cor     = obs_cor,
#     rand_mean   = mean(rand_cors, na.rm=TRUE),
#     rand_sd     = sd(rand_cors,   na.rm=TRUE),
#     perm_pval   = perm_p
#   )
# }
# 
# rank_cor_dt <- rbindlist(rank_cor_list)
# print(rank_cor_dt)
# 
# # ── plot: scatter of discovery rank vs perturb-seq rank per data source ────────
# p <- ggplot(combined[!is.na(rank_pct_perturb)],
#             aes(x = rank_pct_src, y = rank_pct_perturb)) +
#   geom_point(color = '#E41A1C', size = 3, alpha = 0.8) +
#   geom_smooth(method = 'lm', se = TRUE, color = 'grey40',
#               linewidth = 0.8, linetype = 'dashed') +
#   geom_text_repel(aes(label = DPGs), size = 2.5,
#                   max.overlaps = Inf, segment.size = 0.2) +
#   facet_wrap(~ data_source, scales = 'free') +
#   geom_text(data = rank_cor_dt,
#             aes(x = Inf, y = Inf,
#                 label = paste0('r=', round(obs_cor,2),
#                                '\np=', round(perm_pval,3))),
#             inherit.aes = FALSE,
#             hjust = 1.1, vjust = 1.5, size = 3, color = 'grey30') +
#   labs(x        = 'Rank in discovery dataset (%)',
#        y        = 'Rank in CD4 perturb-seq (%)',
#        title    = 'Core gene rank concordance: discovery vs CD4 perturb-seq',
#        subtitle = 'Per data source | Spearman correlation with permutation p-value') +
#   theme_classic(base_size = 11) +
#   theme(plot.title   = element_text(face = 'bold'),
#         strip.text   = element_text(face = 'bold'))
# 
# print(p)
# ggsave('rank_concordance_by_datasource.png', p, width=10, height=5, dpi=150)
# 
# 
# 
# # ── Modification 1: absolute rank in discovery dataset ────────────────────────
# combined[, abs_rank_src := rank_in_src]  # already have this from earlier
# 
# # ── Plot 1: scatter with absolute rank on x-axis ──────────────────────────────
# p1 <- ggplot(combined[!is.na(rank_pct_perturb)],
#              aes(x = abs_rank_src, y = rank_pct_perturb)) +
#   geom_point(color = '#E41A1C', size = 3, alpha = 0.8) +
#   geom_smooth(method = 'lm', se = TRUE, color = 'grey40',
#               linewidth = 0.8, linetype = 'dashed') +
#   geom_text_repel(aes(label = DPGs), size = 2.5,
#                   max.overlaps = Inf, segment.size = 0.2) +
#   facet_wrap(~ data_source, scales = 'free_x') +
#   geom_text(data = rank_cor_dt,
#             aes(x = Inf, y = Inf,
#                 label = paste0('r=', round(obs_cor,2),
#                                '\np=', round(perm_pval,3))),
#             inherit.aes = FALSE,
#             hjust = 1.1, vjust = 1.5, size = 3, color = 'grey30') +
#   scale_x_continuous(breaks = function(x) unique(floor(pretty(x)))) +
#   labs(x        = 'Rank in discovery dataset (1 = most significant)',
#        y        = 'Rank in CD4 perturb-seq (%)',
#        title    = 'Core gene rank concordance: discovery vs CD4 perturb-seq',
#        subtitle = 'Spearman correlation with permutation p-value') +
#   theme_classic(base_size = 11) +
#   theme(plot.title = element_text(face = 'bold'),
#         strip.text = element_text(face = 'bold'))
# 
# print(p1)
# ggsave('rank_concordance_scatter.png', p1, width=10, height=5, dpi=300)
# 
# # ── Plot 2: bar plot of perturb-seq rank % per group ──────────────────────────
# # add random genes — sample once per data source for visualization
# set.seed(42)
# rand_bar_list <- list()
# for (src in data_sources) {
#   core_src   <- combined[data_source == src & !is.na(rank_pct_perturb)]
#   n_core_src <- nrow(core_src)
#   rand_genes <- sample(non_core_perturb, n_core_src, replace = FALSE)
#   rand_dt    <- all_dpg_min[DPG %in% rand_genes,
#                             .(DPGs = DPG, rank_pct_perturb = rank_pct)]
#   rand_dt[, data_source := paste0(src, ' (random)')]
#   rand_bar_list[[src]] <- rand_dt
# }
# rand_bar <- rbindlist(rand_bar_list)
# 
# # combine core + random for bar plot
# bar_plot_dt <- rbind(
#   combined[!is.na(rank_pct_perturb),
#            .(DPGs, data_source, rank_pct_perturb)],
#   rand_bar[, .(DPGs, data_source, rank_pct_perturb)]
# )
# 
# # color: real = red shades per source, random = grey shades
# bar_plot_dt[, group_type := ifelse(grepl('random', data_source),
#                                    'Random', 'Core genes')]
# 
# # mean per group for ordering
# group_means <- bar_plot_dt[, .(mean_pct = mean(rank_pct_perturb)), by = data_source]
# 
# p2 <- ggplot(bar_plot_dt,
#              aes(x = reorder(data_source, rank_pct_perturb,
#                              FUN = mean),
#                  y = rank_pct_perturb,
#                  fill = group_type)) +
#   geom_boxplot(width = 0.5, outlier.size = 0.8, linewidth = 0.6) +
#   geom_jitter(width = 0.1, size = 1.5, alpha = 0.6,
#               aes(color = group_type)) +
#   scale_fill_manual(values  = c('Core genes' = '#E41A1C',
#                                 'Random'      = 'grey70'),
#                     name = NULL) +
#   scale_color_manual(values = c('Core genes' = '#C00000',
#                                 'Random'      = 'grey40'),
#                      guide  = 'none') +
#   scale_y_reverse(labels = function(x) paste0(x, '%')) +
#   labs(x        = 'Group',
#        y        = 'Rank in CD4 perturb-seq (% — lower is better)',
#        title    = 'CD4 perturb-seq rank by discovery dataset group',
#        subtitle = 'Core genes vs matched random genes per data source') +
#   theme_classic(base_size = 11) +
#   theme(plot.title      = element_text(face = 'bold'),
#         legend.position = 'bottom',
#         axis.text.x     = element_text(angle = 15, hjust = 1))
# 
# print(p2)
# ggsave('rank_barplot_by_group.png', p2, width=8, height=5, dpi=150)
# 
# 
# 
# # ── bar plot: proportion of genes in top k rank (absolute) in perturb-seq ─────
# # ── for top k genes in discovery, what % are in top k of perturb-seq ──────────
# top_ks <- c(4, 8, 12, 16)
# 
# bar_count_list <- list()
# for (topk in top_ks) {
#   for (src in data_sources) {
#     
#     # top k core genes in discovery (by rank_in_src)
#     core_src  <- combined[data_source == src & !is.na(rank_pct_perturb)]
#     top_disco <- core_src[rank_in_src <= topk]
#     
#     # how many of those are also in top k% of perturb-seq
#     perturb_topk_pct <- topk / nrow(all_dpg_min) * 100
#     n_in_perturb     <- sum(top_disco$rank_pct_perturb <= perturb_topk_pct)
#     
#     bar_count_list[[paste0(src, '_core_', topk)]] <- data.table(
#       data_source = src,
#       group_type  = 'Core genes',
#       top_k       = topk,
#       n_in        = n_in_perturb,
#       n_total     = nrow(top_disco),
#       proportion  = ifelse(nrow(top_disco) > 0,
#                            n_in_perturb / nrow(top_disco) * 100, NA)
#     )
#     
#     # random: sample same n genes from perturb-seq, check top k% rank
#     rand_props <- replicate(1000, {
#       rand_genes   <- sample(non_core_perturb, topk, replace=FALSE)
#       rand_perturb <- all_dpg_min[DPG %in% rand_genes, rank_pct]
#       mean(rand_perturb <= perturb_topk_pct) * 100
#     })
#     
#     bar_count_list[[paste0(src, '_rand_', topk)]] <- data.table(
#       data_source = paste0(src, ' (random)'),
#       group_type  = 'Random',
#       top_k       = topk,
#       n_in        = NA,
#       n_total     = topk,
#       proportion  = mean(rand_props)
#     )
#   }
# }
# 
# bar_count_dt <- rbindlist(bar_count_list)
# bar_count_dt[, top_k_label := factor(paste0('Top ', top_k, ' in discovery'),
#                                      levels = paste0('Top ', top_ks, ' in discovery'))]
# bar_count_dt[, data_source := factor(data_source,
#                                      levels = c(data_sources,
#                                                 paste0(data_sources, ' (random)')))]
# 
# p <- ggplot(bar_count_dt,
#             aes(x = top_k_label, y = proportion, fill = data_source)) +
#   geom_bar(stat = 'identity', position = position_dodge(0.8), width = 0.7) +
#   geom_text(data = bar_count_dt[group_type == 'Core genes'],
#             aes(label = paste0(n_in, '/', n_total)),
#             position = position_dodge(0.8),
#             vjust = -0.3, size = 2.8) +
#   scale_fill_manual(
#     values = c('DGN'           = '#E41A1C',
#                'eqtl'          = '#377EB8',
#                'DGN (random)'  = '#FCBBA1',
#                'eqtl (random)' = '#BDD7E7'),
#     name = NULL
#   ) +
#   scale_y_continuous(limits = c(0, 115),
#                      breaks = seq(0, 100, 25),
#                      labels = function(x) paste0(x, '%')) +
#   labs(x        = 'Top k genes in discovery dataset',
#        y        = '% also in top k of CD4 perturb-seq',
#        title    = 'Top discovery genes replicate in CD4 perturb-seq',
#        subtitle = 'Core genes (solid) vs random genes (light) per discovery source') +
#   theme_classic(base_size = 11) +
#   theme(plot.title      = element_text(face = 'bold'),
#         legend.position = 'bottom')
# 
# print(p)
# ggsave('rank_topk_discovery_vs_perturb.png', p, width=9, height=5, dpi=150)
# 
# 
# # ── for top k genes in discovery, get their perturb-seq rank % ────────────────
# top_ks <- c(4, 8, 12, 16)
# 
# box_list <- list()
# for (topk in top_ks) {
#   for (src in data_sources) {
#     
#     # top k core genes in discovery
#     core_src  <- combined[data_source == src & !is.na(rank_pct_perturb)]
#     top_disco <- core_src[rank_in_src <= topk]
#     
#     box_list[[paste0(src, '_core_', topk)]] <- data.table(
#       data_source  = src,
#       group_type   = 'Core genes',
#       top_k_label  = paste0('Top ', topk, ' in discovery'),
#       rank_pct_perturb = top_disco$rank_pct_perturb,
#       DPGs         = top_disco$DPGs
#     )
#     
#     # random: sample same n genes, get their perturb-seq rank %
#     rand_genes   <- sample(non_core_perturb, topk, replace=FALSE)
#     rand_perturb <- all_dpg_min[DPG %in% rand_genes,
#                                 .(DPGs = DPG, rank_pct_perturb = rank_pct)]
#     rand_perturb[, data_source := paste0(src, ' (random)')]
#     rand_perturb[, group_type  := 'Random']
#     rand_perturb[, top_k_label := paste0('Top ', topk, ' in discovery')]
#     
#     box_list[[paste0(src, '_rand_', topk)]] <- rand_perturb
#   }
# }
# 
# box_dt <- rbindlist(box_list, fill = TRUE)
# box_dt[, top_k_label := factor(top_k_label,
#                                levels = paste0('Top ', top_ks, ' in discovery'))]
# box_dt[, data_source := factor(data_source,
#                                levels = c(data_sources,
#                                           paste0(data_sources, ' (random)')))]
# 
# p <- ggplot(box_dt, aes(x = top_k_label, y = rank_pct_perturb,
#                         fill = data_source)) +
#   geom_boxplot(width = 0.6, outlier.size = 0.8, linewidth = 0.6,
#                position = position_dodge(0.8)) +
#   geom_jitter(aes(color = data_source),
#               position = position_jitterdodge(jitter.width = 0.1,
#                                               dodge.width  = 0.8),
#               size = 1.8, alpha = 0.8) +
#   geom_text_repel(
#     data     = box_dt[group_type == 'Core genes' & !is.na(DPGs)],
#     aes(label = DPGs, color = data_source),
#     position = position_jitterdodge(jitter.width=0.1, dodge.width=0.8),
#     size = 2.2, max.overlaps = Inf, segment.size = 0.2, show.legend = FALSE
#   ) +
#   scale_fill_manual(
#     values = c('DGN'           = '#E41A1C',
#                'eqtl'          = '#377EB8',
#                'DGN (random)'  = '#FCBBA1',
#                'eqtl (random)' = '#BDD7E7'),
#     name = NULL
#   ) +
#   scale_color_manual(
#     values = c('DGN'           = '#C00000',
#                'eqtl'          = '#1A5276',
#                'DGN (random)'  = '#E8A898',
#                'eqtl (random)' = '#85B0C8'),
#     guide = 'none'
#   ) +
#   scale_y_reverse(labels = function(x) paste0(x, '%')) +
#   labs(x        = 'Top k genes in discovery dataset',
#        y        = 'Rank in CD4 perturb-seq (% — lower is better)',
#        title    = 'Perturb-seq rank of top discovery genes',
#        subtitle = 'Core genes (solid) vs random (light) | lower % = higher rank') +
#   theme_classic(base_size = 11) +
#   theme(plot.title      = element_text(face = 'bold'),
#         legend.position = 'bottom')
# 
# p <- ggplot(box_dt, aes(x = top_k_label, y = rank_pct_perturb,
#                         fill = data_source)) +
#   geom_boxplot(width = 0.6, outlier.size = 0.8, linewidth = 0.6,
#                position = position_dodge(0.8)) +
#   geom_jitter(aes(color = data_source),
#               position = position_jitterdodge(jitter.width = 0.1,
#                                               dodge.width  = 0.8),
#               size = 1.8, alpha = 0.8) +
#   scale_fill_manual(
#     values = c('DGN'           = '#E41A1C',
#                'eqtl'          = '#377EB8',
#                'DGN (random)'  = '#FCBBA1',
#                'eqtl (random)' = '#BDD7E7'),
#     name = NULL
#   ) +
#   scale_color_manual(
#     values = c('DGN'           = '#C00000',
#                'eqtl'          = '#1A5276',
#                'DGN (random)'  = '#E8A898',
#                'eqtl (random)' = '#85B0C8'),
#     guide = 'none'
#   ) +
#   scale_y_reverse(labels = function(x) paste0(x, '%')) +
#   labs(x        = 'Top k genes in discovery dataset',
#        y        = 'Rank in CD4 perturb-seq (% — lower is better)',
#        title    = 'Perturb-seq rank of top discovery genes',
#        subtitle = 'Core genes (solid) vs random (light) | lower % = higher rank') +
#   theme_classic(base_size = 11) +
#   theme(plot.title      = element_text(face = 'bold'),
#         legend.position = 'bottom')
# 
# print(p)
# ggsave('rank_topk_boxplot_perturb.png', p, width=10, height=6, dpi=300)
# 
# 
# 
