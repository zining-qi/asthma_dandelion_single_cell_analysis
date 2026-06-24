# Asthma Dandelion Single Cell Analysis

## Overview

This repository accompanies the study:

*Leveraging trans-gene regulation prioritizes central genes and pathways in asthma.*


This repository focuses on **downstream characterization and validation of disease-proximal genes (DPGs)** identified using the DANDELION framework.

The analysis consists of two main components:

1. **Single-cell analysis**
   We leverage human lung cell atlas (HLCA; Sikkema et al., 2023 Nat Med) single-cell RNA-seq data to characterize the **cell-type-specific activity and expression patterns of DPGs**, including:

   * Quantification of gene activity using AUCell
   * Identification of relevant cell types for DPG
   * Visualization of expression patterns across cell populations (DotPlot and UMAP)

2. **Replication analysis in CD4+ Perturb-seq data**
   We validate DPG prioritization by performing replication analysis in **primary CD4+ Perturb-seq datasets** (Zhu et al., 2025, bioRxiv).



The workflow consists of two main components:

1. **Single-cell analysis** to evaluate cell-type-specific expression and activity of DPGs in lung tissue
2. **Replication analysis** of DPGs using CD4+ Perturb-seq data

---

## Repository Structure

```
asthma_dandelion_single_cell_analysis/
├── asthma_dandelion_replication/
│   ├── asthma_replication_analysis.R
│   └── dandelion_cd4_asthma.R
├── asthma_single_cell_analysis/
│   ├── asthma_dpg_lung_atlas_aucell_calculation.R
│   ├── asthma_dpg_lung_atlas_aucell_celltype_nomination.R
│   ├── asthma_dpg_lung_atlas_sc_expr_dotplot.R
│   └── asthma_dpg_lung_atlas_sc_expr_umap.R
└── README.md
```

---

## 1. Single Cell Analysis (Human Lung Cell Atlas)

**Directory:** `asthma_single_cell_analysis/`

This module analyzes the activity and expression of DPGs in **HLCA**.

### Files

* **`asthma_dpg_lung_atlas_aucell_calculation.R`**
  Computes **AUCell scores** across cell type based on DPG to quantify gene activity.

* **`asthma_dpg_lung_atlas_aucell_celltype_nomination.R`**
  Performs **cell-type nomination** based on AUCell scores to identify relevant cell populations for DPG.

* **`asthma_dpg_lung_atlas_sc_expr_dotplot.R`**
  Generates **dot plots** showing expression patterns of DPGs across cell types.

* **`asthma_dpg_lung_atlas_sc_expr_umap.R`**
  Produces **UMAP visualizations** of single-cell data with DPG expression and cell-type annotations.

---

## 2. Dandelion Replication Analysis

**Directory:** `asthma_dandelion_replication/`

This module performs replication of DPGs using **primary CD4+ Perturb-seq data**.

### Files

* **`dandelion_cd4_asthma.R`**
  Implements the **Dandelion pipeline** to identify disease-proximal genes in asthma.

* **`asthma_replication_analysis.R`**
  Performs **post-Dandelion replication analysis**, including evaluation of identified DPGs in perturbation data.

---


---

## Notes

* Input data paths are currently hard-coded and should be updated for portability.
* Large single-cell datasets may require substantial memory.
* Gene naming consistency (e.g., SYMBOL vs ENSG) is critical across all analyses.

---

## Citation

If you use this code or find it helpful for your research, please cite:

Salamone, I. M.†, Tian, P.†, Qi, Z., Zhao, J., Zhang, L., Tan, Q., Li, J., Michael, A. N., Thornburg, A. G., Sakabe, N. J., Weber, Z. T., Minogue, M., Chen, B., Ciszewski, C., He, X., Shah, H., Vercelli, D., Ober, C., Lin, H., Liu, Z.‡, Nóbrega, M. A.‡, & Liu, X.‡ (under review).
**Leveraging trans-gene regulation prioritizes central genes and pathways in asthma.**

† These authors contributed equally to this work.
‡ Corresponding authors.

