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
DimPlot(MPE.combined, reduction = "umap")
