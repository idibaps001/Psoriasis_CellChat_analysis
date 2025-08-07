#Load the required libraries
library(CellChat)
library(patchwork)
options(future.globals.maxSize= 8912896000)

#Load cellchat objects
setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/01_z1_wt_untreated/cd45negANDcd45pos/cellchat_v2")
cellchat.Z1WTUNT <- readRDS("cellchat.Z1WTUNT.rds")

setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/02_z1_wt_imq/cd45negANDcd45pos/cellchat_v2")
cellchat.Z1WTIMQ <- readRDS("cellchat.Z1WTIMQ.rds")

#Set a directory to save figures
setwd("/Users/lucia/Documents/23.psoriasis/08.CellchatComparison/01.Z1WTIMQvsZ1WTUNT/cellchat_v2/")
opar = par(no.readonly = TRUE)

#Merge two cellchat object togather
object.list <- list(Z1WTUNT = cellchat.Z1WTUNT, Z1WTIMQ = cellchat.Z1WTIMQ)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
#An object of class CellChat created from a merged object with multiple datasets 
#1245 signaling genes.
#13530 cells. 
#CellChat analysis of single cell RNA-seq data! 

### Part I: Predict general principles of cell-cell communication
#Compare the total number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),size.text = 16, color.use = c("#5d6663","#ba963f" ))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",size.text = 16, color.use = c("#5d6663","#ba963f" ))
p = gg1 + gg2
ggsave(p, file = "barplot1.tiff", height = 5, width = 7)

saveRDS(cellchat, file = "cellchat_comparisonAnalysis_Z1WTUNT_vs_Z1WTIMQ.rds")

#Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations

tiff("Differential number of interactions.tiff", width = 1500, height = 1500)
netVisual_diffInteraction(cellchat, weight.scale = T, title.name = "Differential number of interactions",vertex.label.cex = 3, edge.width.max = 24,comparison = c(1, 2),margin = 0.5, 
                          color.use =  c( "Dendritic epidermal T cells" =  "#8DD3C7",
                                          "Mono/Macrophages" = "#FFFFB3",
                                          "T/NK cells" = "#BEBADA",
                                          "DCs"   = "#FB8072",
                                          "Neutrophils"  = "#80B1D3",
                                          "Epidermal cells" = "#FDB462",
                                          "Stromal cells" = "#B3DE69",
                                          "Endothelial cells" = "#FCCDE5")   )
dev.off() 

tiff("Differential interaction strength.tiff", width = 1500, height = 1500)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", title.name = "Differential interaction strength",vertex.label.cex = 3, edge.width.max = 24,comparison = c(1, 2),margin = 0.5, 
                          color.use =  c( "Dendritic epidermal T cells" =  "#8DD3C7",
                                          "Mono/Macrophages" = "#FFFFB3",
                                          "T/NK cells" = "#BEBADA",
                                          "DCs"   = "#FB8072",
                                          "Neutrophils"  = "#80B1D3",
                                          "Epidermal cells" = "#FDB462",
                                          "Stromal cells" = "#B3DE69",
                                          "Endothelial cells" = "#FCCDE5")   )
dev.off() 

#heatmap
tiff("Differential number of interactions_heatmap.tiff", width = 5, height = 5, units = "in", res = 300)
library(RColorBrewer)
display.brewer.all()  
brewer.pal(8, "Set3")  
# [1] "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5"
netVisual_heatmap(cellchat, comparison = c(1,2),font.size = 10, font.size.title = 12, 
                   color.use =  c( "Dendritic epidermal T cells" =  "#8DD3C7",
                                  "Mono/Macrophages" = "#FFFFB3",
                                  "T/NK cells" = "#BEBADA",
                                  "DCs"   = "#FB8072",
                                  "Neutrophils"  = "#80B1D3",
                                  "Epidermal cells" = "#FDB462",
                                  "Stromal cells" = "#B3DE69",
                                  "Endothelial cells" = "#FCCDE5")   )

#[1] "Dendritic epidermal T cells" "Mono/Macrophages"            "T/NK cells"                  "DCs"                        
#[5] "Neutrophils"                 "Epidermal cells"             "Stromal cells"               "Endothelial cells"        
dev.off() #红色代表第2组的细胞互作数目更多，互作更频繁 （Z1WTIMQ组是第2组）

#> Do heatmap based on a merged object
tiff("Differential interaction strength_heatmap.tiff", width = 5, height = 5, units = "in", res = 300)
 netVisual_heatmap(cellchat, measure = "weight", comparison = c(1,2),font.size = 10, font.size.title = 12, 
                   color.use =  c( "Dendritic epidermal T cells" =  "#8DD3C7",
                                   "Mono/Macrophages" = "#FFFFB3",
                                   "T/NK cells" = "#BEBADA",
                                   "DCs"   = "#FB8072",
                                   "Neutrophils"  = "#80B1D3",
                                   "Epidermal cells" = "#FDB462",
                                   "Stromal cells" = "#B3DE69",
                                   "Endothelial cells" = "#FCCDE5")   )
 dev.off() 


##Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
tiff("Relative information flow_Z1WTIMQvsZ1WTUNT.tiff", width = 9, height = 13, units = "in", res = 300)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(2, 1), font.size = 8)
dev.off() 
#Black colors mean equally important, one color like red means the enriched pathways in the first condition, and the other color like blue means the enriched pathways in the second condition. From the bar size, you can also make a judgement if the difference is large.
#Bar size represents the information flow, which is computed by summarizing all the communication probability in that signaling pathway.

setwd("/Users/lucia/Documents/23.psoriasis/08.CellchatComparison/01.Z1WTIMQvsZ1WTUNT/cellchat_v3/")
tiff("Relative information flow_Z1WTIMQvsZ1WTUNT_SelectedSources.tiff", width = 5, height = 7, units = "in", res = 300)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(2, 1), font.size = 8, sources.use = "Mono/Macrophages")
dev.off() 

setwd("/Users/lucia/Documents/23.psoriasis/08.CellchatComparison/01.Z1WTIMQvsZ1WTUNT/cellchat_v3/")
tiff("Relative information flow_Z1WTIMQvsZ1WTUNT_SelectedSourcesTargets.tiff", width = 4, height = 5, units = "in", res = 300)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(2, 1), font.size = 8, sources.use = "Mono/Macrophages", targets.use = "Epidermal cells")
dev.off() 

par(opar)
