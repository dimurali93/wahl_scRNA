source("~/dimurali/wahl_scRNA/scripts/library.R")
dir.create("/home/divya/dimurali/wahl_scRNA/seurat")
setwd("/home/divya/dimurali/wahl_scRNA/seurat")

# read the snRNA assay from the 10x raw matrix(the filtering is a little different from 10X, some cells may not be there)

counts.adult <- Read10X("/home/divya/dimurali/wahl/cellranger_results/K18_663_SinglePosEpi_P1/outs/filtered_feature_bc_matrix")
counts.d3to5tp <- Read10X("/home/divya/dimurali/wahl/cellranger_results/scrna_d3to5tp/outs/filtered_feature_bc_matrix")
counts.w8tp <- Read10X("/home/divya/dimurali/wahl/cellranger_results/scrna_w8tp/outs/filtered_feature_bc_matrix/")

# create seurat object
adult <- CreateSeuratObject(counts.adult, project = "Adult")
d3to5tp <- CreateSeuratObject(counts.d3to5tp, project = "D3to5")
w8tp <- CreateSeuratObject(counts.w8tp, project = "W8")

#make a list of all seurat objects and name them
scRNA.list <- list(d3to5tp, w8tp, adult) 
names(scRNA.list) <- c( "d3to5tp", "w8tp", "adult") 

#markers list
markers_extendedlist <- c("Krt8","Krt18","Sox11","Sox4","Hspa1b","Hspa1a","Col1a1","Dcn","Ptprc","Csf1r","Upk3b","Kdr","Upk3a","Cdh19")
markers <- c("Trp63","Krt14","Krt5","Kit","Ehf","Elf5","Mki67","Foxa1","Esr1","Pgr","Prlr")

#preporcess data such that as naming Idents, FindVariableFeatures, ScaleData, RunPCA, cumsum, RunUMAP, FindNeighbors,FindClusters
#nfeatures=2000
source("~/dimurali/wahl_scRNA/scripts/functions/processNewSeurat.R")
scRNA.list <- lapply(scRNA.list, processNewSeurat)

#PREQC feature, ncount, nfeature, passedQC plots
source("~/dimurali/wahl_scRNA/scripts//functions/processPlots.R")
path="/home/divya/dimurali/wahl_scRNA/seurat/preQC/"
dir.create(path)
scRNA.list[[1]] <- processPlots(scRNA.list[[1]], "d3to5_preQC",path)
scRNA.list[[2]] <- processPlots(scRNA.list[[2]], "w8_preQC",path)
scRNA.list[[3]] <- processPlots(scRNA.list[[3]], "adult_preQC",path)


## filter dead cells out of all clusters and run scaling and normalisation again
source("~/dimurali/wahl_scRNA/scripts/functions/filterDeadCells.R")
scRNA.list <- lapply(scRNA.list, filterDeadCells)
scRNA.list <- lapply(scRNA.list, processNewSeurat)
# Length before filtering
# [1] 37958
# Length after filtering
# [1] 34115
# Number of dead cells = 3843
# Length before filtering
# [1] 5569
# Length after filtering
# [1] 3783
# Number of dead cells = 1786
# Length before filtering
# [1] 6386
# Length after filtering
# [1] 5393
# Number of dead cells = 993


source("~/dimurali/wahl_scRNA/scripts/functions/genPostFilteredFeaturePlots.R")
dir.create(path)
path="/home/divya/dimurali/wahl_scRNA/seurat/genPostFilteredFeaturePlots/"
scRNA.list[[1]] <- genPostFilteredFeaturePlots(scRNA.list[[1]], "d3to5_filtered", markers, path)
scRNA.list[[2]] <- genPostFilteredFeaturePlots(scRNA.list[[2]], "w8_filtered",markers, path)
scRNA.list[[3]] <- genPostFilteredFeaturePlots(scRNA.list[[3]], "adult",markers, path)
# extended markers1 containing others as well
scRNA.list[[1]] <- genPostFilteredFeaturePlots(scRNA.list[[1]], "d3to5_filtered_ext",markers_extendedlist, path)
scRNA.list[[2]] <- genPostFilteredFeaturePlots(scRNA.list[[2]], "w8_filtered_ext",markers_extendedlist, path)
scRNA.list[[3]] <- genPostFilteredFeaturePlots(scRNA.list[[3]], "adult_ext",markers_extendedlist, path)

# Display the filtered plots
d3to5_preQC_clusters=DimPlot(scRNA.list[[1]], label = T)+ plot_annotation(title = 'D3to5 clusters')
w8_preQC_clusters=DimPlot(scRNA.list[[2]], label = T) + plot_annotation(title = 'W8 clusters')
Adult_preQC_clusters=DimPlot(scRNA.list[[3]], label = T) + plot_annotation(title = 'Adult')

AllSamples_dotplot <- ggarrange(d3to5_preQC_clusters, w8_preQC_clusters, Adult_preQC_clusters,
                        nrow = 1, ncol = 3)
ggsave("AllSamples_dotplot.png")
dev.off()

