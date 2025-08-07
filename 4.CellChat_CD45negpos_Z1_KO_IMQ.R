#CellChat 
library(Seurat)
library(tidyverse)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 8912896000)

setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/03_z1_ko_imq/cd45pos/")
load("CD45POS_Z1_KO_IMQ.rds")
CD45POS_Z1_KO_IMQ <- Z1_KO_IMQ

setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/03_z1_ko_imq/cd45negANDcd45pos/")
load("CD45NEG_Z1_KO_IMQ.rds")
CD45NEG_Z1_KO_IMQ <- Z1_KO_IMQ

#Merge normalized seurat objects
#https://satijalab.org/seurat/articles/merge_vignette.html
CD45POS_Z1_KO_IMQ <- NormalizeData(CD45POS_Z1_KO_IMQ)
CD45NEG_Z1_KO_IMQ <- NormalizeData(CD45NEG_Z1_KO_IMQ)
Z1_KO_IMQ <- merge(CD45POS_Z1_KO_IMQ, y = CD45NEG_Z1_KO_IMQ,add.cell.ids = c("cd45pos", "cd45neg"), project = "psoriasis", merge.data = TRUE)

length(Cells(CD45POS_Z1_KO_IMQ)) #[1] 4568
length(Cells(CD45NEG_Z1_KO_IMQ)) #[1] 4766
length(Cells(Z1_KO_IMQ)) #[1] 9334


###################################### Modify the names of the three cells to be consistent with the manuscript version.
Idents(Z1_KO_IMQ) <- "celltypes_2"
levels(Z1_KO_IMQ)
#"Dendritic epidermal T cells", "M1-like Macrophages", "Foxp3+ Tregs", "Vγ6 T cells",  "Vγ4 T cells",  "Slamf9+ Macrophages", "Conventional CD4+ T cells", "Langerhans cells",  "Selenop+ M2-like Macrophages", "DoubleNegative T cells", "Resident M2-like Macrophages", "Monocytes", "ILC2", "Neutrophils", "MHCII+Ccr2+ Macrophages", "Cytotoxic CD8+ T cells", "Mature Monocytes" , "Proliferative CD4+ T cells", "Proliferative LC",  "pDC", "Granular keratinocytes", "Stromal cells cycling", "Basal keratinocytes", "Supra-Basal keratinocytes", "Universal fibroblasts", "Basal keratinocytes cycling", "Endothelial cells", "Reticular fibroblasts", "Inner ear fibroblasts", "Hair dermal sheath", "Cornified keratinocytes", "Dermal fibrocytes"    

Z1_KO_IMQ@meta.data$celltypes_4 <- case_when(  
  Z1_KO_IMQ@meta.data$celltypes_2 == "Selenop+ M2-like Macrophages" ~ "Selenop+ Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_2 == "Resident M2-like Macrophages" ~ "Resident Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_2 == "M1-like Macrophages"  ~ "Thbs1+ Macrophages",
  TRUE ~ Z1_KO_IMQ@meta.data$celltypes_2
)

#check
table( Z1_KO_IMQ@meta.data$celltypes_2,Z1_KO_IMQ@meta.data$celltypes_4    )
table( Z1_KO_IMQ@meta.data$celltypes_2    )
table( Z1_KO_IMQ@meta.data$celltypes_4    ) 

#select cells
Idents(Z1_KO_IMQ) <- "celltypes_4"
levels(Z1_KO_IMQ)

Z1_KO_IMQ <- subset(Z1_KO_IMQ, idents = c("Dendritic epidermal T cells", "Thbs1+ Macrophages", "Foxp3+ Tregs", "Vγ6 T cells",  "Vγ4 T cells",  "Slamf9+ Macrophages", 
                                          "Conventional CD4+ T cells", "Langerhans cells",  "Selenop+ Macrophages", "DoubleNegative T cells", "Resident Macrophages",
                                          "Monocytes", "ILC2", "Neutrophils", "MHCII+Ccr2+ Macrophages", "Cytotoxic CD8+ T cells", "Mature Monocytes" , "Proliferative CD4+ T cells", 
                                          "Proliferative LC", "pDC", "Granular keratinocytes", "Stromal cells cycling", "Basal keratinocytes", "Supra-Basal keratinocytes", "Universal fibroblasts",
                                          "Basal keratinocytes cycling", "Endothelial cells", "Reticular fibroblasts", "Inner ear fibroblasts", "Hair dermal sheath", "Cornified keratinocytes", "Dermal fibrocytes"    
))
length(Cells(Z1_KO_IMQ)) #[1] 7649 cells are selected

# chage the level of naming system
Z1_KO_IMQ@meta.data$celltypes_5 <- case_when(  
  Z1_KO_IMQ@meta.data$celltypes_4 == "Thbs1+ Macrophages" ~ "Mono/Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Foxp3+ Tregs" ~ "T/NK cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Vγ6 T cells"  ~ "T/NK cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Vγ4 T cells"  ~ "T/NK cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Slamf9+ Macrophages" ~ "Mono/Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Conventional CD4+ T cells" ~ "T/NK cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Langerhans cells"  ~ "DCs",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Selenop+ Macrophages" ~ "Mono/Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_4 == "DoubleNegative T cells" ~ "T/NK cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Resident Macrophages"  ~ "Mono/Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Monocytes" ~ "Mono/Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_4 == "ILC2" ~ "T/NK cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "MHCII+Ccr2+ Macrophages"  ~ "Mono/Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Cytotoxic CD8+ T cells" ~ "T/NK cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Mature Monocytes" ~ "Mono/Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Proliferative CD4+ T cells"  ~ "Mono/Macrophages",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Proliferative LC" ~ "DCs",
  Z1_KO_IMQ@meta.data$celltypes_4 == "pDC" ~ "DCs",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Granular keratinocytes"  ~ "Epidermal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Stromal cells cycling" ~ "Stromal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Basal keratinocytes" ~ "Epidermal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Supra-Basal keratinocytes"  ~ "Epidermal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Universal fibroblasts" ~ "Stromal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Basal keratinocytes cycling" ~ "Epidermal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Reticular fibroblasts"  ~ "Stromal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Inner ear fibroblasts" ~ "Stromal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Hair dermal sheath" ~ "Stromal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Cornified keratinocytes"  ~ "Epidermal cells",
  Z1_KO_IMQ@meta.data$celltypes_4 == "Dermal fibrocytes" ~ "Stromal cells",
  
  TRUE ~ Z1_KO_IMQ@meta.data$celltypes_4
)

table( Z1_KO_IMQ@meta.data$celltypes_5, Z1_KO_IMQ@meta.data$celltypes_4   ) %>% view()
table( Z1_KO_IMQ@meta.data$celltypes_5  ) 

Idents(Z1_KO_IMQ) <- "celltypes_5"
levels(Z1_KO_IMQ)
#[1] "Mono/Macrophages"            "Dendritic epidermal T cells" "T/NK cells"                  "DCs"                        
#[5] "Neutrophils"                 "Stromal cells"               "Endothelial cells"           "Epidermal cells"   


#CellChat tutorial https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html
data.input <- GetAssayData(Z1_KO_IMQ, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Z1_KO_IMQ)
samples <- as.factor( Z1_KO_IMQ@meta.data$conditions )
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

#add metadata and group info
#cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
#cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
#check
unique(meta$group)
levels(cellchat@idents) 
#[1] "Mono/Macrophages"            "Dendritic epidermal T cells" "T/NK cells"                  "DCs"                        
#[5] "Neutrophils"                 "Stromal cells"               "Endothelial cells"           "Epidermal cells"      

table(cellchat@idents)
# Mono/Macrophages            Dendritic epidermal T cells                  T/NK cells                         DCs 
#           2145                         101                                    847                           129 
# Neutrophils               Stromal cells             Endothelial cells             Epidermal cells 
#    531                        3073                         179                         644 

#CellChat tutorial https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
#showDatabaseCategory(CellChatDB)
#Show the structure of the database
#dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis 
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = 2) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

#generate a table of cell communication analysis results
df.net <- subsetCommunication(cellchat)# a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors
df.netP <- subsetCommunication(cellchat,slot.name = "netP") # the inferred communications at the level of signaling pathways

dim(df.net)
#[1] 1635    11
dim(df.netP)
#[1] 823     5

setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/03_z1_ko_imq/cd45negANDcd45pos/cellchat_v2/")
write.csv2(df.net, file = "CD45NEGPOS-Z1_KO_IMQ-CellChat-net_v2.csv")
write.csv2(df.netP, file = "CD45NEGPOS-Z1_KO_IMQ-CellChat-netP_v2.csv")

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  
  }


#Part III: Visualization of cell-cell communication network
cellchat@netP$pathways

pathways.show <- c( "TNF"  )
vertex.receiver = seq(1,2)  # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.TNF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TNF[1,] # show one ligand-receptor pair
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show)

netVisual_bubble(cellchat, sources.use = 20, targets.use = c(3:8), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 7, targets.use = c(13:20), remove.isolate = FALSE)

#Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "TNF")

#Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
#Compute and visualize the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 20, height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 20, height = 20)
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht3 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("TNF","IL4","IL6"),pattern = "outgoing",width = 20, height = 20)
ht4 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("TNF","IL4","IL6"),pattern = "incoming",width = 20, height = 20)
ht3 + ht4

#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
#Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 2#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 2.
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, height = 15)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

#Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming") #Cophenetic values begin to drop when the number of incoming patterns is 4.
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, height = 15)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")


###################################################################################################
#                                                                                                 #
#                     assign a new name to each CellChat dataset                                  #
#                                                                                                 #
#                             cellchat.Z1KOIMQ <- cellchat                                        #
#                                                                                                 #
###################################################################################################
setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/03_z1_ko_imq/cd45negANDcd45pos/cellchat_v2/")
saveRDS(cellchat, file = "cellchat.Z1KOIMQ.rds")
