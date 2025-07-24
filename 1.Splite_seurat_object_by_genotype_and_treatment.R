#The Seurat object `cd45pos` was partitioned into six datasets based on genotype and treatment. 
#These datasets are designated as follows: Z1_WT_UNTREATED, Z1_WT_IMQ, Z1_KO_IMQ, Z2_WT_UNTREATED, Z2_WT_IMQ, and Z2_KO_IMQ.
setwd("/Users/lucia/Documents/23.psoriasis/01.Psoriasis_2023-5-25")
library(Seurat)
library(tidyverse)
load("post12_cd45pos_alldata_final_clustering_annotation.rds")
Idents(post12_cd45pos_alldata_final_clustering_annotation) <- "condition"
levels(post12_cd45pos_alldata_final_clustering_annotation)
#[1] "Zeb1-WT_IMQ"  "Zeb1-KO_IMQ"  "Zeb2-WT_IMQ"  "Zeb2-KO_IMQ"  "Zeb1-WT_None" "Zeb2-WT_None"

Z1_WT_UNTREATED <- subset(post12_cd45pos_alldata_final_clustering_annotation, idents = c( "Zeb1-WT_None") )
Z1_WT_IMQ <- subset(post12_cd45pos_alldata_final_clustering_annotation, idents = c( "Zeb1-WT_IMQ") )
Z1_KO_IMQ <- subset(post12_cd45pos_alldata_final_clustering_annotation, idents = c( "Zeb1-KO_IMQ") )

setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat")
save( Z1_WT_UNTREATED   , file = "01_z1_wt_untreated/cd45pos/CD45POS_Z1_WT_UNTREATED.rds")
save( Z1_WT_IMQ   , file = "02_z1_wt_imq/cd45pos/CD45POS_Z1_WT_IMQ.rds")
save( Z1_KO_IMQ   , file = "03_z1_ko_imq/cd45pos/CD45POS_Z1_KO_IMQ.rds")

length(Cells(Z1_WT_UNTREATED ))
#[1] 2996
length(Cells(Z1_WT_IMQ ))
#[1] 4720
length(Cells(Z1_KO_IMQ ))
#[1] 4568


#The Seurat object `cd45eg` was partitioned into six datasets based on genotype and treatment. 
#These datasets are designated as follows: Z1_WT_UNTREATED, Z1_WT_IMQ, Z1_KO_IMQ, Z2_WT_UNTREATED, Z2_WT_IMQ, and Z2_KO_IMQ.
setwd("/Users/lucia/Documents/23.psoriasis/01.Psoriasis_2023-5-25")
library(Seurat)
library(tidyverse)
load("post12_cd45neg_alldata_final_clustering_annotation.rds")
Idents(post12_cd45neg_alldata_final_clustering_annotation) <- "condition"
levels(post12_cd45neg_alldata_final_clustering_annotation)
#[1] "Zeb1-WT_IMQ"  "Zeb1-KO_IMQ"  "Zeb2-WT_IMQ"  "Zeb2-KO_IMQ"  "Zeb1-WT_None" "Zeb2-WT_None"

Z1_WT_UNTREATED <- subset(post12_cd45neg_alldata_final_clustering_annotation, idents = c( "Zeb1-WT_None") )
Z1_WT_IMQ <- subset(post12_cd45neg_alldata_final_clustering_annotation, idents = c( "Zeb1-WT_IMQ") )
Z1_KO_IMQ <- subset(post12_cd45neg_alldata_final_clustering_annotation, idents = c( "Zeb1-KO_IMQ") )

setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat")
save( Z1_WT_UNTREATED   , file = "01_z1_wt_untreated/cd45negANDcd45pos/CD45NEG_Z1_WT_UNTREATED.rds")
save( Z1_WT_IMQ   , file = "02_z1_wt_imq/cd45negANDcd45pos/CD45NEG_Z1_WT_IMQ.rds")
save( Z1_KO_IMQ   , file = "03_z1_ko_imq/cd45negANDcd45pos/CD45NEG_Z1_KO_IMQ.rds")

length(Cells(Z1_WT_UNTREATED ))
#[1] 3861
length(Cells(Z1_WT_IMQ ))
#[1] 5235
length(Cells(Z1_KO_IMQ ))
#[1] 4766
