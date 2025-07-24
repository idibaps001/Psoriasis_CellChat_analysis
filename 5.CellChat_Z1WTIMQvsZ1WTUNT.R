#Load the required libraries
library(CellChat)
library(patchwork)
#Set a directory to save figures
setwd("/Users/lucia/Documents/23.psoriasis/08.CellchatComparison/01.Z1WTIMQvsZ1WTUNT")
opar = par(no.readonly = TRUE)

#Merge two cellchat object togather
object.list <- list(Z1WTUNT = cellchat.Z1WTUNT, Z1WTIMQ = cellchat.Z1WTIMQ)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
#An object of class CellChat created from a merged object with multiple datasets 
#887 signaling genes.
#13530 cells. 
#CellChat analysis of single cell RNA-seq data! 

################################## Suppl Fig 2F ##################################
### Part I: Predict general principles of cell-cell communication
#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

saveRDS(cellchat, file = "cellchat_comparisonAnalysis_Z1WTUNT_vs_Z1WTIMQ.rds")

#Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#heatmap
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object = object.list[[i]])
  }
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

####Part II: Identify the conserved and context-specific signaling pathways
###Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
##Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional",do.parallel = FALSE)
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2

################################## Suppl Fig 2G ##################################
##Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
#Black colors mean equally important, one color like red means the enriched pathways in the first condition, and the other color like blue means the enriched pathways in the second condition. From the bar size, you can also make a judgement if the difference is large.
#Bar size represents the information flow, which is computed by summarizing all the communication probability in that signaling pathway.


par(opar)
