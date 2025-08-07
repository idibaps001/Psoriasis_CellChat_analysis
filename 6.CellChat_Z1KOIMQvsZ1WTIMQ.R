#Load the required libraries
library(CellChat)
library(patchwork)
options(future.globals.maxSize= 8912896000)

#Load cellchat objects
setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/02_z1_wt_imq/cd45negANDcd45pos/cellchat_v2")
cellchat.Z1WTIMQ <- readRDS("cellchat.Z1WTIMQ.rds")

setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/03_z1_ko_imq/cd45negANDcd45pos/cellchat_v2")
cellchat.Z1KOIMQ <- readRDS("cellchat.Z1KOIMQ.rds")

#Set a directory to save figures
setwd("/Users/lucia/Documents/23.psoriasis/08.CellchatComparison/02.Z1KOIMQvsZ1WTIMQ/cellchat_v2/")
opar = par(no.readonly = TRUE)

#Merge two cellchat object togather
object.list <- list(Z1WTIMQ = cellchat.Z1WTIMQ, Z1KOIMQ = cellchat.Z1KOIMQ)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
#An object of class CellChat created from a merged object with multiple datasets 
#1245 signaling genes.
#15659 cells. 
#CellChat analysis of single cell RNA-seq data! 

### Part I: Predict general principles of cell-cell communication
#Compare the total number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),size.text = 16, color.use = c("#ba963f","#95594c" ))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",size.text = 16, color.use = c("#ba963f","#95594c" ))
p = gg1 + gg2
ggsave(p, file = "barplot1.tiff", height = 5, width = 7)
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_Z1KOIMQ_vs_Z1WTIMQ.rds")

#Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations
tiff("Differential number of interactions.tiff", width = 1500, height = 1500)
netVisual_diffInteraction(cellchat, weight.scale = T, title.name = "Differential number of interactions",vertex.label.cex = 3, edge.width.max = 24,comparison = c(1, 2),margin = 0.3, 
                          color.use =  c( "#FFFFB3", "#8DD3C7" , "#FB8072", "#80B1D3","#BEBADA", "#FDB462" , "#FCCDE5" ,"#B3DE69")   )
dev.off() 

tiff("Differential interaction strength.tiff", width = 1500, height = 1500)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", title.name = "Differential interaction strength",vertex.label.cex = 3, edge.width.max = 24,comparison = c(1, 2),margin = 0.3, 
                          color.use =  c( "#FFFFB3", "#8DD3C7" , "#FB8072", "#80B1D3","#BEBADA", "#FDB462" , "#FCCDE5" ,"#B3DE69")   )
dev.off() 

#heatmap
tiff("Differential number of interactions_heatmap.tiff", width = 5, height = 5, units = "in", res = 300)
netVisual_heatmap(cellchat, comparison = c(1,2),font.size = 10, font.size.title = 12, 
                  color.use =  c( "#FFFFB3", "#8DD3C7" , "#FB8072", "#80B1D3","#BEBADA", "#FDB462" , "#FCCDE5" ,"#B3DE69")   )
dev.off() 

#> Do heatmap based on a merged object
tiff("Differential interaction strength_heatmap.tiff", width = 5, height = 5, units = "in", res = 300)
netVisual_heatmap(cellchat, measure = "weight", comparison = c(1,2),font.size = 10, font.size.title = 12, 
                  color.use =  c( "#FFFFB3", "#8DD3C7" , "#FB8072", "#80B1D3","#BEBADA", "#FDB462" , "#FCCDE5" ,"#B3DE69")   )
dev.off() 

##Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
tiff("Relative information flow_Z1WTIMQvsZ1WTUNT.tiff", width = 9, height = 13, units = "in", res = 300)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(2, 1), font.size = 8)
dev.off() 
#Black colors mean equally important, one color like red means the enriched pathways in the first condition, and the other color like blue means the enriched pathways in the second condition. From the bar size, you can also make a judgement if the difference is large.
#Bar size represents the information flow, which is computed by summarizing all the communication probability in that signaling pathway.

setwd("/Users/lucia/Documents/23.psoriasis/08.CellchatComparison/02.Z1KOIMQvsZ1WTIMQ/cellchat_v3/")
tiff("Relative information flow_Z1WTIMQvsZ1WTUNT_SelectedSources.tiff", width = 5, height = 7, units = "in", res = 300)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(2, 1), font.size = 8, sources.use = "Mono/Macrophages")
dev.off() 


par(opar)
