suppressMessages(library(scRepertoire))
library(Seurat) 

# Data load 
MPE.combined_label2 <- readRDS("G:/SingleCell/singlecellfile/MPE-CD4_MPE.combined_label2.rds")
SKOM002 <- read.csv("/SingleCell/KOM002/KOM002_vdj2/outs/filtered_contig_annotations.csv")
SLK208 <- read.csv("/SingleCell/LK208/LK208_vdj_outs/filtered_contig_annotations.csv")
SLK006 <- read.csv("/SingleCell/LK006/outs/filtered_contig_annotations.csv")

contig_list_MPE <- list(SKOM002,SLK208,SLK006)
head(contig_list_MPE[[2]])

combined_MPE <- combineTCR(contig_list_MPE,samples = c("KOM002","LK208","LK006"),ID = c("T","T","T"), cells = "T-AB")

seurat0.8 <- combineExpression(combined_MPE, MPE.combined_label2, 
                               cloneCall="gene", proportion = FALSE, 
                               cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

DimPlot(seurat0.8, group.by = "cloneType",pt.size = 1.5) 

#<Save the object at this point so that it can easily be loaded back>
saveRDS(seurat0.8, file = "G:/SingleCell/singlecellfile/MPE-CD4_seurat0.8.rds")
