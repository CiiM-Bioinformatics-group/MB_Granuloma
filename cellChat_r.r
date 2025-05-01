library(anndata)
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(zellkonverter)
library(SingleCellExperiment)
library(future)
# Load the libraries
library(gridExtra)
library(patchwork)
library(tidyverse)
# Increase the limit to 2 GB or any size you find suitable
options(future.globals.maxSize = 10*2 * 1024^3)  # 2 GB * 10

# reticulate::use_python("/vol/projects/jxun/.R_env/bin/python", required=T) 
# ad <- read_h5ad(("/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Anndata/AllDataHarmolized_01.h5ad"))

sce <- zellkonverter::readH5AD(file = "/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Anndata/G_niches.h5ad")
# retrieve all the available assays within sce objec

str(sce)

# extract a cell meta data
assay(sce, "logcounts") <- assay(sce, "X")
meta <- as.data.frame(SingleCellExperiment::colData(sce)) #
# cellchat <- createCellChat(object = sce, group.by = "Spots")

cellchat <- createCellChat(object = sce, group.by = "Spots")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)


# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692
# print(as.numeric(execution.time, units = "secs"))
#> [1] 13.20763
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
# print(as.numeric(execution.time, units = "secs"))

rm(list=ls())
load("/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Anndata/G_niches.Rdata")
ls()

# head(cellchat@net)|


cellchat@idents <- recode(cellchat@idents,
                      "Layer2 (Fibro & Plasma Cells)" = "Outer Layer",
                      "Layer1 (Plasma Cells & Macro & T & Fibro)" = "Middle Layer",
                      "Layer0 (Macro & Mono & NK & T)" = "Inner Layer",
                     "Necrotizing" = "Necrotizing",
                     "T & B" = "T & B"
                    )
meta <- cellchat@meta
meta$Spots <- recode(meta$Spots,
                      "Layer2 (Fibro & Plasma Cells)" = "Outer Layer",
                      "Layer1 (Plasma Cells & Macro & T & Fibro)" = "Middle Layer",
                      "Layer0 (Macro & Mono & NK & T)" = "Inner Layer",
                     "Necrotizing" = "Necrotizing",
                     "T & B" = "T & B"
                    )
cellchat@meta <- meta

# Define your custom labels for the rows and columns
custom_labels <- c("Necrotizing", "Inner Layer", "Middle Layer", "Outer Layer", "T & B")
# Update row and column names for cellchat@net$count
rownames(cellchat@net$count) <- custom_labels
colnames(cellchat@net$count) <- custom_labels
# Update row and column names for cellchat@net$weight
rownames(cellchat@net$weight) <- custom_labels
colnames(cellchat@net$weight) <- custom_labels

# Update the row and column names of the prob array in netP
dimnames(cellchat@netP$prob)[[1]] <- custom_labels
dimnames(cellchat@netP$prob)[[2]] <- custom_labels

# Update the clusters column in var.features
cellchat@var.features$features.info$clusters <- recode(cellchat@var.features$features.info$clusters,
                                                        "Layer2 (Fibro & Plasma Cells)" = "Outer Layer",
                                                        "Layer1 (Plasma Cells & Macro & T & Fibro)" = "Middle Layer",
                                                        "Layer0 (Macro & Mono & NK & T)" = "Inner Layer",
                                                        "Necrotizing" = "Necrotizing",
                                                        "T & B" = "T & B")
cellchat <- setIdent(cellchat, ident.use = "Spots")

# Update the labels in cellchat@idents to match the new cluster names
levels(cellchat@idents) <- custom_labels

rgb_to_hex <- function(r, g, b) {
  rgb(r, g, b, maxColorValue = 1)
}
# Define colors in R
group.colors <- c(
  "Necrotizing" = rgb_to_hex(1.0, 0.596078431372549, 0.5882352941176471), # Uncomment if needed
  "Inner Layer" = rgb_to_hex(0.5803921568627451, 0.403921568627451, 0.7411764705882353),
  "Middle Layer" = rgb_to_hex(0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
  # "AT & Plasma Cells" = rgb_to_hex(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), # Uncomment if needed
  "Outer Layer" = rgb_to_hex(0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
  "T & B" = rgb_to_hex(1.0, 0.7333333333333333, 0.47058823529411764)
)
group.colors

options(repr.plot.width = 8, repr.plot.height = 6)
groupSize <- as.numeric(table(cellchat@idents))
head(cellchat@meta)
# Define custom colors for each group
#group.colors <- c("Necrotizing"="", "Inner Layer" = "red", "Middle Layer" = "blue", "Outer Layer" = "green", "T & B" = "purple")
# par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", vertex.label.cex = F, color.use = group.colors
                 )


netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, vertex.label.cex = 1.5,color.use = group.colors
                )

# save(cellchat, file="/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Anndata/G_niches.Rdata")
# load("/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Anndata/G_niches.Rdata")
ls()

options(repr.plot.width =36, repr.plot.height = 7)
mat <- cellchat@net$weight
par(mfrow = c(1,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, edge.weight.max = max(mat),  vertex.label.cex = 4, color.use = group.colors)#, #title.name = rownames(mat)[i])
}

mat <- cellchat@net$count
par(mfrow = c(1,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), vertex.label.cex = 4,color.use = group.colors)#, #title.name = rownames(mat)[i])
}


options(repr.plot.width = 6, repr.plot.height = 5.5)
pathways.show <- c("MIF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
# vertex.receiver = seq(1,4) # a numeric vector. 
# netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.label.cex = 1.5, color.use = group.colors)

# Chord diagram
options(repr.plot.width = 8, repr.plot.height = 8)
pathways.show =  c("MIF") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = group.colors)

options(repr.plot.width = 8, repr.plot.height = 22)
pathways.show <- c("TGFb", 
                   "CXCL", 
                   "SPP1", 
                   "IL16", "CD99",
                   "CCL", "MIF", "ICAM", 
                   "COLLAGEN", 
                   "COMPLEMENT",
                   "APP", "ApoE", "PECAM1", "BAFF", "APRIL", "CD96") 
# pathways.show <- c("COLLAGEN")
netAnalysis_contribution(cellchat, signaling = pathways.show)

options(repr.plot.width = 8, repr.plot.height = 8)
par(mfrow = c(2,4), xpd=TRUE)
pairLR <- c("APP_CD74", "SPP1_CD44", "MIF_CD74_CD44", "CXCL12_CXCR4", "CD99_PILRA", "ICAM1_ITGAX_ITGB2", "IL16_CD4", "CCL5_CCR1")
for (i in 1:length(pairLR)){
    # pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
    LR.show <- pairLR[i] # show one ligand-receptor pair
    # Hierarchy plot
    netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show)
    #> [[1]]
    # Circle plot
    netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
}

# Load the library
library(cowplot)

pairLR

options(repr.plot.width = 6, repr.plot.height = 6)
plots <- list()
# Assuming you have already set up the signaling pathways to display
pathways.show <- c("COLLAGEN")
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show)$interaction_name
#pairLR <- c("APP_CD74", "SPP1_CD44", "MIF_CD74_CD44", "CXCL12_CXCR4", "CD99_PILRA", "ICAM1_ITGAX_ITGB2", "IL16_CD4", "CCL5_CCR1")
# Create individual plots and store them in the list
for (i in 1:length(pairLR)) {
    LR.show <- pairLR[i]  # Show one ligand-receptor pair
    # Circle plot
    p2 <- netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, 
                               layout = "circle", vertex.label.cex = 1.6, color.use = group.colors)
    plots[i] <- p2
}

length(plots)

combined_plot <- plot_grid(
    plotlist = plots,  # List of plots
    nrow = 2,          # Number of rows
    ncol = 2,           # Number of columns
    align = "hv",         # Align both horizontally and vertically
    axis = "tblr",        # Ensure axes are aligned properly
    # rel_heights = c(2, 2),# Relative heights for the two rows (adjust if needed)
    # rel_widths = rep(1, 4),# Relative widths for the four columns
    labels = NULL,        # No additional labels to avoid clutter
    label_size = 10,      # Optional: label size for the subplots
    label_y = 1.1,
    scale = 0.95  
)


options(repr.plot.width = 8, repr.plot.height = 8)
combined_plot

# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = c(1:5), targets.use = c(1:5), 
                 signaling = c("MIF", "COMPLEMENT", "ICAM",  "CXCL"), 
                 remove.isolate = TRUE)
bubble_plot_flipped
#> Comparing communications on a single object

options(repr.plot.width = 7, repr.plot.height = 24) 
# Call netVisual_bubble() with custom colors
bubble_plot <- netVisual_bubble(
  cellchat,
  sources.use = c(1:5),
  targets.use = c(1:5),
  signaling = c("MIF", "SPP1", "ICAM", "CD99", "CXCL", 
                #"GALECTIN", "PECAM1", 
                #"ApoE", "CCL", 
                "IL16", "LAIR1", "TGFb", "BAFF", "APRIL", 
                #"CD45", "NOTCH", "CD6", 
                #"CD40",  "CD23",  
                "TNF", "IL2", 
                #"TWEAK", "CD80", "CD86", 
                "TRAIL", "MHC-I",  #"CD39", 
                "CD96" 
                #"FLT3"
               ),
  remove.isolate = TRUE
  
)

# Modify the plot to swap the x and y axes
bubble_plot_flipped <- bubble_plot + coord_flip() + 
theme(
  text = element_text(size = 8)  # Reduce the font size
)

options(repr.plot.width = 14, repr.plot.height = 5) 
bubble_plot_flipped

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
options(repr.plot.width = 10, repr.plot.height = 10) 
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5), targets.use = c(1,2,3,4,5), lab.cex = 0.45, 
                     signaling = c("MIF", "CXCL"),small.gap=TRUE,
                     legend.pos.y = 10)

options(repr.plot.width = 16, repr.plot.height = 12) 
# Extract all communication data for a specific signaling pathway or the entire dataset
communication_data <- subsetCommunication(cellchat, slot.name = "net")
communication_data <- dplyr::filter(communication_data, annotation=="ECM-Receptor")
# Sort the data by the strength of interactions or another metric, then select the top N interactions
# Replace "contribution" with the appropriate column that indicates the interaction strength or contribution
top_communication_data <- communication_data[order(communication_data$prob, decreasing = TRUE), ]
top_communication_data <- head(top_communication_data, 50)  # Replace 10 with the desired number of top pairs

# Extract the ligand-receptor pairs from the top communication data
top_LR_pairs <- unique(paste(top_communication_data$ligand, top_communication_data$receptor, sep = "_"))
top_LR_pairs_df <- data.frame(interaction_name = top_LR_pairs, stringsAsFactors = FALSE)

table(communication_data$pathway_name)
head(communication_data)

# dplyr::filter(communication_data, ligand=="COL6A1")

# dplyr::filter(communication_data, ligand=="COL1A1")

# Use netVisual_bubble() to visualize the selected top L-R pairs
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5), targets.use = c(1,2,3,4,5), lab.cex = 1.5, 
                     pairLR.use = top_LR_pairs_df, title="ECM", 
                     legend.pos.y = 20, legend.pos.x = 0, color.use = group.colors)


# Extract all communication data for a specific signaling pathway or the entire dataset
communication_data <- subsetCommunication(cellchat, slot.name = "net")
table(communication_data$annotation)
communication_data <- dplyr::filter(communication_data, annotation=="Cell-Cell Contact")
# Sort the data by the strength of interactions or another metric, then select the top N interactions
# Replace "contribution" with the appropriate column that indicates the interaction strength or contribution
top_communication_data <- communication_data[order(communication_data$prob, decreasing = TRUE), ]
top_communication_data <- head(top_communication_data, 50)  # Replace 10 with the desired number of top pairs

# Extract the ligand-receptor pairs from the top communication data
top_LR_pairs <- unique(paste(top_communication_data$ligand, top_communication_data$receptor, sep = "_"))
top_LR_pairs_df <- data.frame(interaction_name = top_LR_pairs, stringsAsFactors = FALSE)
# Use netVisual_bubble() to visualize the selected top L-R pairs
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5), targets.use = c(1,2,3,4,5), lab.cex = 1.2, 
                     pairLR.use = top_LR_pairs_df, #title="Cell-Cell Contact",
                     legend.pos.y = 20, legend.pos.x = 0, color.use = group.colors)

# Extract all communication data for a specific signaling pathway or the entire dataset
communication_data <- subsetCommunication(cellchat, slot.name = "net")
table(communication_data$annotation)
communication_data <- dplyr::filter(communication_data, annotation=="Secreted Signaling")
# Sort the data by the strength of interactions or another metric, then select the top N interactions
# Replace "contribution" with the appropriate column that indicates the interaction strength or contribution
top_communication_data <- communication_data[order(communication_data$prob, decreasing = TRUE), ]
top_communication_data <- head(top_communication_data, 50)  # Replace 10 with the desired number of top pairs

# Extract the ligand-receptor pairs from the top communication data
top_LR_pairs <- unique(paste(top_communication_data$ligand, top_communication_data$receptor, sep = "_"))
top_LR_pairs_df <- data.frame(interaction_name = top_LR_pairs, stringsAsFactors = FALSE)
# Use netVisual_bubble() to visualize the selected top L-R pairs
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5), targets.use = c(1,2,3,4,5), lab.cex = 1.2, 
                     pairLR.use = top_LR_pairs_df, #title="Secreted Signaling",
                     legend.pos.y = 20, legend.pos.x = 0,color.use = group.colors)

options(repr.plot.width = 5, repr.plot.height = 5) 
plots <- list()
for (i in 1:1){
    pathways.show = c( "COLLAGEN", "ICAM","CXCL", "MIF", "COMPLEMENT", "SPP1")[i]
    cat(i)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
    # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
    p = netAnalysis_signalingRole_network(cellchat, 
                                          signaling = pathways.show, 
                                          width = 5, height = 3, font.size = 12, 
                                          color.use = group.colors)
    print(p)
    plots[i] <- p
}

plots[1]

library(cowplot)
combined_plot <- plot_grid(
    plotlist = plots,  # List of plots
    nrow = 3,          # Number of rows
    ncol = 2,           # Number of columns
    align = "hv",         # Align both horizontally and vertically
    axis = "tblr",        # Ensure axes are aligned properly
    # rel_heights = c(2, 2),# Relative heights for the two rows (adjust if needed)
    # rel_widths = rep(1, 4),# Relative widths for the four columns
    labels = NULL,        # No additional labels to avoid clutter
    label_size = 10,      # Optional: label size for the subplots
    label_y = 1.1,
    scale = 0.95  
)
combined_plot

pathways.show =  c("SPP1") 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 12)

pathways.show =  c("ICAM") 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 12)

options(repr.plot.width = 7, repr.plot.height = 7) 
pathways.show =  c("COMPLEMENT") 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 12)



# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "MIF", "SPP1", "COLLAGEN"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1
gg2

options(repr.plot.width = 10, repr.plot.height = 30) 
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 6, height = 50, font.size = 16, font.size.title = 20, color.use = group.colors)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 6, height = 50, font.size = 16, font.size.title = 20, color.use = group.colors)


library(NMF)
library(ggalluvial)

options(repr.plot.width = 7, repr.plot.height = 7) 
selectK(cellchat, pattern = "outgoing")

options(repr.plot.width = 12, repr.plot.height = 18) 
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = c("outgoing", "incoming"), 
                                          k = nPatterns, color.use = group.colors, font.size = 12,
                                          width = 6,
                                          height = 6)

p <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, color.use = group.colors)

# river plot
options(repr.plot.width = 9, repr.plot.height = 6.5) 
netAnalysis_river(cellchat, pattern = "incoming", color.use = group.colors) #incoming
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
options(repr.plot.width = 16, repr.plot.height = 2.5) 
netAnalysis_dot(cellchat, pattern = "incoming", color.use = group.colors)

netVisual_bubble(cellchat, 
                 sources.use = c(2), 
                 targets.use = c(3,4,5), remove.isolate = TRUE, 
                 signaling = c("COLLAGEN","ICAM","CXCL","MIF","COMPLEMENT","SPP1"))

reticulate::use_python("/home/jxun/.cache/R/basilisk/1.10.2/zellkonverter/1.8.0/zellkonverterAnnDataEnv-0.8.0/bin/python", required=T) 
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

options(repr.plot.width = 6, repr.plot.height = 5) 
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds", font.size = 15, color.use = group.colors)
#> Do heatmap based on a single object

options(repr.plot.width = 4, repr.plot.height = 3) 
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds", font.size = 10, color.use = group.colors)




