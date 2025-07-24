#CellChat 
library(Seurat)
library(tidyverse)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/01_z1_wt_untreated/cd45pos/")
load("CD45POS_Z1_WT_UNTREATED.rds")
CD45POS_Z1_WT_UNTREATED <- Z1_WT_UNTREATED

setwd("/Users/lucia/Documents/23.psoriasis/07.cellchat/01_z1_wt_untreated/cd45negANDcd45pos/")
load("CD45NEG_Z1_WT_UNTREATED.rds")
CD45NEG_Z1_WT_UNTREATED <- Z1_WT_UNTREATED

#Merge normalized seurat objects
#https://satijalab.org/seurat/articles/merge_vignette.html
CD45POS_Z1_WT_UNTREATED <- NormalizeData(CD45POS_Z1_WT_UNTREATED)
CD45NEG_Z1_WT_UNTREATED <- NormalizeData(CD45NEG_Z1_WT_UNTREATED)
Z1_WT_UNTREATED <- merge(CD45POS_Z1_WT_UNTREATED, y = CD45NEG_Z1_WT_UNTREATED,add.cell.ids = c("cd45pos", "cd45neg"), project = "psoriasis", merge.data = TRUE)

length(Cells(CD45POS_Z1_WT_UNTREATED)) #[1] 2996
length(Cells(CD45NEG_Z1_WT_UNTREATED)) #[1] 3861
length(Cells(Z1_WT_UNTREATED)) #[1] 6857

Idents(Z1_WT_UNTREATED) <- "celltypes_2"
levels(Z1_WT_UNTREATED)
#"Dendritic epidermal T cells", "M1-like Macrophages", "Foxp3+ Tregs", "Vγ6 T cells",  "Vγ4 T cells",  "Slamf9+ Macrophages", "Conventional CD4+ T cells", "Langerhans cells",  "Selenop+ M2-like Macrophages", "DoubleNegative T cells", "Resident M2-like Macrophages", "Monocytes", "ILC2", "Neutrophils", "MHCII+Ccr2+ Macrophages", "Cytotoxic CD8+ T cells", "Mature Monocytes" , "Proliferative CD4+ T cells", "Proliferative LC",  "pDC", "Granular keratinocytes", "Stromal cells cycling", "Basal keratinocytes", "Supra-Basal keratinocytes", "Universal fibroblasts", "Basal keratinocytes cycling", "Endothelial cells", "Reticular fibroblasts", "Inner ear fibroblasts", "Hair dermal sheath", "Cornified keratinocytes", "Dermal fibrocytes"    
Z1_WT_UNTREATED <- subset(Z1_WT_UNTREATED, idents = c("Dendritic epidermal T cells", "M1-like Macrophages", "Foxp3+ Tregs", "Vγ6 T cells",  "Vγ4 T cells",  "Slamf9+ Macrophages", "Conventional CD4+ T cells", "Langerhans cells",  "Selenop+ M2-like Macrophages", "DoubleNegative T cells", "Resident M2-like Macrophages", "Monocytes", "ILC2", "Neutrophils", "MHCII+Ccr2+ Macrophages", "Cytotoxic CD8+ T cells", "Mature Monocytes" , "Proliferative CD4+ T cells", "Proliferative LC", "pDC", "Granular keratinocytes", "Stromal cells cycling", "Basal keratinocytes", "Supra-Basal keratinocytes", "Universal fibroblasts", "Basal keratinocytes cycling", "Endothelial cells", "Reticular fibroblasts", "Inner ear fibroblasts", "Hair dermal sheath", "Cornified keratinocytes", "Dermal fibrocytes"    
))


#CellChat Tutorial https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html
data.input <- GetAssayData(Z1_WT_UNTREATED, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(Z1_WT_UNTREATED)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

#add metadata and set default cell identity
#cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
#cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
#check grouping info
unique(meta$group)
levels(cellchat@idents) 
#[1] "Dendritic epidermal T cells"  "M1-like Macrophages"          "Foxp3+ Tregs"                 "Vγ6 T cells"                  "Vγ4 T cells"                  "Slamf9+ Macrophages"          "Conventional CD4+ T cells"    "Langerhans cells"             "Selenop+ M2-like Macrophages"
#[10] "DoubleNegative T cells"       "Resident M2-like Macrophages" "Monocytes"                    "ILC2"                         "Neutrophils"                  "MHCII+Ccr2+ Macrophages"      "Cytotoxic CD8+ T cells"       "Mature Monocytes"             "Proliferative CD4+ T cells"  
#[19] "Proliferative LC"             "pDC"                          "Granular keratinocytes"       "Stromal cells cycling"        "Basal keratinocytes"          "Supra-Basal keratinocytes"    "Universal fibroblasts"        "Basal keratinocytes cycling"  "Endothelial cells"           
#[28] "Reticular fibroblasts"        "Inner ear fibroblasts"        "Hair dermal sheath"           "Cornified keratinocytes"      "Dermal fibrocytes"           
table(cellchat@idents)
#Dendritic epidermal T cells    M1-like Macrophages       Foxp3+ Tregs                  Vγ6 T cells                  Vγ4 T cells          Slamf9+ Macrophages    Conventional CD4+ T cells 
#1305                           29                           27                          127                          117                           28                           48 
#Langerhans cells Selenop+ M2-like Macrophages       DoubleNegative T cells Resident M2-like Macrophages        Monocytes                         ILC2                  Neutrophils 
#43                          107                          129                           53                           11                           89                           12 
#MHCII+Ccr2+ Macrophages       Cytotoxic CD8+ T cells    Mature Monocytes   Proliferative CD4+ T cells             Proliferative LC                pDC 
#74                            9                           13                           12                            2                             1 
#Granular keratinocytes        Stromal cells cycling     Basal keratinocytes    Supra-Basal keratinocytes        Universal fibroblasts  Basal keratinocytes cycling       Endothelial cells 
#285                           41                          358                           87                          784                          244                          184 
#Reticular fibroblasts        Inner ear fibroblasts        Hair dermal sheath      Cornified keratinocytes            Dermal fibrocytes 
#781                          226                          111                          165                           18 
#CellChat Tutorial https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html


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
future::plan("multicore", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

#generate cell-cell communication results
df.net <- subsetCommunication(cellchat)# a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors
df.netP <- subsetCommunication(cellchat,slot.name = "netP") # the inferred communications at the level of signaling pathways

dim(df.net)
#[1] 13904    11
dim(df.netP)
#[1] 6482    5

write.csv2(df.net, file = "CD45NEGPOS-Z1_WT_UNTREATED-CellChat-net.csv")
write.csv2(df.netP, file = "CD45NEGPOS-Z1_WT_UNTREATED-CellChat-netP.csv")


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
#[1] "MIF"         "CD45"        "TNF"         "APP"         "ICAM"        "CCL"         "THY1"        "FN1"         "CDH1"        "MHC-II"      "THBS"        "GDF"         "CD48"       
#[14] "GALECTIN"    "CXCL"        "CD96"        "CD137"       "CDH"         "BST2"        "SPP1"        "CD86"        "NECTIN"      "IGF"         "ITGAL-ITGB2" "TGFb"        "LAMININ"    
#[27] "SELPLG"      "CD40"        "TWEAK"       "NKG2D"       "SEMA4"       "VISFATIN"    "MHC-I"       "APRIL"       "COMPLEMENT"  "ALCAM"       "CD6"         "CD52"        "SEMA6"      
#[40] "CLEC"        "OSM"         "ANNEXIN"     "IL1"         "GAS"         "ICOS"        "IL2"         "CD80"        "CD226"       "PVR"         "IL10"        "CD39"        "LAIR1"      
#[53] "JAM"         "LIGHT"       "L1CAM"       "IL4"         "IL6"         "OCLN"        "RANKL"       "PD-L1"       "MPZ"         "PDL2"       

pathways.show <- c( "TNF"  )
vertex.receiver = seq(1,2)  # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.TNF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TNF[1,] # show one ligand-receptor pair
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "hierarchy")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "spatial")

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
#   To facilitate intergroup comparisons, each CellChat dataset was assigned a new name.          #
#                                                                                                 #
#                             cellchat.Z1WTUNT <- cellchat                                        #
#                                                                                                 #
###################################################################################################
