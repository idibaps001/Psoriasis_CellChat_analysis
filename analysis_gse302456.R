# ------------------------------------------------------------------------------
# Title: Macrophage Analysis in Psoriasis (GSE302456)
# Description: Differential expression and subclustering analysis of macrophages
#              from scRNA-seq data of psoriatic and healthy skin
# Original dataset by Chen Q, Wei L, Liu Y, Li Y, Jia X, Li Y, Du X
# Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302456
# ------------------------------------------------------------------------------
# Clear workspace
rm(list=ls())

# Load required libraries
# R version 4.5.0 (2025-04-11)
library(Seurat) # 5.3.0
library(patchwork) # 1.3.0
library(ggplot2) # 3.5.2

# Load Seurat object from RDS file
setwd("/Users/albertoatencia/Desktop/psoriasis/GSE302456/")
rds_path <- "data/GSE302456_GCBS0428822.rds"
seurat_gse302456 <- readRDS(rds_path)

# ------------------------------------------------------------------------------
# 1. Visualize original clustering by celltype
DimPlot(seurat_gse302456, group.by = "celltype", label=T) + ggtitle("Original Clustering by Cell Type")

# ------------------------------------------------------------------------------
# 2. Subset macrophage cells
macrophage_cells <- WhichCells(seurat_gse302456, expression = celltype == "Macrophages")
macrophage_obj <- subset(seurat_gse302456, cells = macrophage_cells)
# Macrophages count by group
print(table(Group = macrophage_obj@meta.data$group))

# Differential expression: PS vs HS
Idents(macrophage_obj) <- "group"    # Set cell identities based on group metadata - HS(Healthy Scalp), PS (Psoriatic Scalp), PB (Psoriatic Body)
deg_PS_vs_HS <- FindMarkers(macrophage_obj, ident.1 = "PS", ident.2 = "HS", logfc.threshold = 0.25, min.pct = 0.1)

# Visualize genes of interest
DotPlot(macrophage_obj, 
        features = c("CYBA", "TXN", "HMOX1", "FABP5", "FTL", "FTH1", "PRDX1", "MRC1", "GAS5"), 
        scale = TRUE) +
  coord_flip() +
  scale_color_distiller(palette = "RdYlBu", direction = -1) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +  
  theme_minimal() +
  labs(title = "Macrophages", x = NULL, y = NULL) +
  theme(
    axis.text = element_text(size = 12),
    axis.text.x = element_text(face = "italic"),  
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.background = element_blank(),
    panel.grid.major = element_blank()
  )

# ------------------------------------------------------------------------------