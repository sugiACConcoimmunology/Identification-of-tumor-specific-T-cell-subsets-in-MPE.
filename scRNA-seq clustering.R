library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratData)

dfKOM002 <- Read10X(data.dir = "filtered_feature_bc_matrix")
dfLK208 <- Read10X(data.dir = "filtered_feature_bc_matrix")
dfLK006 <- Read10X(data.dir = "filtered_feature_bc_matrix")

# Calculate TRBV abundances on the raw counts before creating a Seurat object
KOM002_TRB.index <- grep(pattern = "^TRB", x = rownames(dfKOM002), value = F) 
LK208_TRB.index <- grep(pattern = "^TRB", x = rownames(dfLK208), value = F) 
LK006_TRB.index <- grep(pattern = "^TRB", x = rownames(dfLK006), value = F) 
# Select row indices and not TRBV names 
# Remove TRBV from count.data
dfKOM002 <- dfKOM002[-KOM002_TRB.index, ]
dfLK208 <- dfLK208[-LK208_TRB.index, ]
dfLK006 <- dfLK006[-LK006_TRB.index, ]
# Calculate TRAV abundances on the raw counts before creating a Seurat object
KOM002_TRA.index <- grep(pattern = "^TRA", x = rownames(dfKOM002), value = F) 
LK208_TRA.index <- grep(pattern = "^TRA", x = rownames(dfLK208), value = F) 
LK006_TRA.index <- grep(pattern = "^TRA", x = rownames(dfLK006), value = F) 
# Select row indices and not TRAV names 
# Remove TRAV from count.data
dfKOM002 <- dfKOM002[-KOM002_TRA.index, ]
dfLK208 <- dfLK208[-LK208_TRA.index, ]
dfLK006 <- dfLK006[-LK006_TRA.index, ]


# Create Seurat object
ddKOM002 <- CreateSeuratObject(counts = dfKOM002,
                               project = "KOM002",
                               min.cells = 3, 
                               min.features = 200, 
                               row.names = 1)
ddLK006 <- CreateSeuratObject(counts = dfLK006,
                                project = "LK006",
                                min.cells = 3, 
                                min.features = 200, 
                                row.names = 1)
ddLK208 <- CreateSeuratObject(counts = dfLK208,
                              project = "LK208",
                              min.cells = 3, 
                              min.features = 200, 
                              row.names = 1)

ddKOM002 <- subset(ddKOM002, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10) 
ddLK208 <- subset(ddLK208, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
ddLK006 <- subset(ddLK006, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Merge cases and remove batch effects
ddMPE <- merge(ddKOM002,y=c(ddLK208,ddLK006),add.cell.ids=c("KOM002_T","LK208_T","LK006_T"))
head(colnames(ddMPE))
tail(colnames(ddMPE))
unique(sapply(X = strsplit(colnames(ddMPE), split = "_"), FUN = "[", 1))
table(ddMPE$orig.ident)
ddMPE

ddMPE.list <- SplitObject(ddMPE,split.by="orig.ident")
ddMPE.list <- lapply(X = ddMPE.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
ddMPE.features <- SelectIntegrationFeatures(object.list = ddMPE.list)
MPE.anchors <- FindIntegrationAnchors(object.list = ddMPE.list, anchor.features = ddMPE.features,k.anchor = 20)

MPE.combined <- IntegrateData(anchorset = MPE.anchors)
DefaultAssay(MPE.combined) <- "integrated"

MPE.combined <- ScaleData(MPE.combined, verbose = FALSE)
MPE.combined <- RunPCA(MPE.combined, npcs = 30, verbose = FALSE)
MPE.combined <- RunUMAP(MPE.combined, reduction = "pca", dims = 1:30)
MPE.combined <- FindNeighbors(MPE.combined, reduction = "pca", dims = 1:30)
MPE.combined <- FindClusters(MPE.combined, resolution = 0.4)
DimPlot(MPE.combined, reduction = "umap" , label = TRUE,pt.size = 1.5, label.size = 7) 
DimPlot(MPE.combined, reduction = "umap" , label.size = 3,split.by = "orig.ident")

cluster.all.markers <- FindAllMarkers(MPE.combined,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Top10_MPE <- cluster.all.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
DoHeatmap(MPE.combined,features = Top10_MPE$gene)

# Averaged feature expression by identity class
MPE.averages <- AverageExpression(MPE.combined)
head(MPE.averages[["RNA"]][, 1:5])
orig.levels <- levels(MPE.combined)
Idents(MPE.combined) <- gsub(pattern = " ", replacement = "_", x = Idents(MPE.combined))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(MPE.combined) <- orig.levels
MPE.averages <- AverageExpression(MPE.combined, return.seurat = TRUE)
MPE.averages

# Naming representative marker
Tex <- c("PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","ENTPD1","ITGAE","ICOS","TNFRSF9")
Tact <- c("CD38","CD69","ICOS","TNFRSF9","HLA-DRA","CD40LG")
Tef <- c("GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY","IFNG","FASLG","TNF","IL17A","CD69","FGFBP2","KLRD1","S1PR1")
Tn <- c("SELL","CCR7","IL7R","CD28","FAS","ITGAL","S1PR1","LEF1","SATB1")
Tnkgd  <- c("MKI67","TCF7","TOX","TRGV9","TRDV2","KLRB1","KLRC3","XCL2","NKG7","NCR1")

DoHeatmap(MPE.averages, features = c(Tn,Tnkgd),slot = "scale.data" ,size = 4,draw.lines = FALSE) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
DoHeatmap(MPE.averages, features = c(Tex4,Tact,Tef,"MKI67"),slot = "scale.data" ,size = 4,draw.lines = FALSE) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))

# Cluster labeling
new.cluster.ids <- c("C4-CX3CR1","C2-IL7R","C1-LEF1","C5-PDCD1","C3-RPL13","C6-ZNF683","C8-LMNA","C10-S100A9","C7-MKI67","C9-IFIT2")
new.cluster.ids2 <- c("C4","C2","C1","C5","C3","C6","C8","C10","C7","C9")
names(new.cluster.ids) <- levels(MPE.combined)
MPE.combined_label <- RenameIdents(MPE.combined , new.cluster.ids)
names(new.cluster.ids2) <- levels(MPE.combined)
MPE.combined_label2 <- RenameIdents(MPE.combined , new.cluster.ids2)

#<Save the object at this point so that it can easily be loaded back>
saveRDS(MPE.combined, file = "G:/SingleCell/singlecellfile/MPE_MPE.combined.rds")
saveRDS(MPE.combined_label, file = "G:/SingleCell/singlecellfile/MPE_MPE.combined_label.rds")
saveRDS(MPE.combined_label2, file = "G:/SingleCell/singlecellfile/MPE_MPE.combined_label2.rds")

