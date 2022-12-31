setwd("/home/divya/dimurali/wahl")
library("dplyr")
library("Seurat")
library("gbutils")
library("ggplot2")
processSeurat <-function(object,
           nfeatures = 2000,
           dims = 1:30,
           resolution = 0.5,
           normalize = T) {
    if (normalize) {
      x <-
        NormalizeData(object) %>% 
        FindVariableFeatures(nfeatures = nfeatures) %>% 
        ScaleData() %>% 
        RunPCA() %>% 
        RunUMAP(dims =dims) %>% 
        FindNeighbors(dims = dims) %>% 
        FindClusters(resolution = resolution)
    } else{
      x <-
        FindVariableFeatures(object) %>% 
        ScaleData(nfeatures = nfeatures) %>% 
        RunPCA() %>% RunUMAP(dims =dims) %>% 
        FindNeighbors(dims = dims) %>% 
        FindClusters(resolution = resolution)
    }
  }
# import cellranger output and create seurat object
d3to5tp <- Read10X(data.dir = "/home/divya/dimurali/wahl/results/scrna_d3to5tp/outs/filtered_feature_bc_matrix")
# convert into seurat object
d3to5tp <- CreateSeuratObject(counts = d3to5tp,
                                     project = "d3to5tp", min.cells = 3, min.features = 200)
# %mt cells
d3to5tp[["percent.mt"]] <- PercentageFeatureSet(d3to5tp, pattern = "^mt-")
tail(rownames(d3to5tp))
vlnplot_d3to5tp=VlnPlot(d3to5tp, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
ggsave("./plot_preQC/vlnplot_d3to5tp.png",width = 8,height = 8)
# name the orig.ident to use when merging all time points
d3to5tp$orig.ident <- "d3to5tp"
# basic processing and normalisation
d3to5tp<- processSeurat(d3to5tp)
# umap1
dimplot_d3to5tp=DimPlot(d3to5tp, label = T)+labs(color = "Clusters")
ggsave("./plot_preQC/dimplot_d3to5tp.png",width = 8,height = 8)
# markers for annotation
markers <-
  c(
    "Krt5",
    "Krt14",
    "Trp63",
    "Acta2",
    "Krt8",
    "Elf5",
    "Ehf",
    "Kit",
    "Esr1",
    "Pgr",
    "Foxa1",
    "Mki67"
  )
featureplot_d3to5tp_markers=FeaturePlot(d3to5tp,
            features = markers,
            ncol = 4,
            label = T)# LOOKS GOOD BASAL TO LUM TRANSITION
ggsave("./plot_preQC/featureplot_d3to5tp.png",width = 8,height = 8)
FeaturePlot(
  d3to5tp_seurat,
  features = c("percent.mt", "Mki67", "nCount_RNA"),
  ncol = 3,
  label = T
)#CLUSTER 6 IS HIGH with mt-genes
d3to5tp_seurat <- FindClusters(d3to5tp_seurat, resolution = 0.8)
dimplot_d3to5tp=DimPlot(d3to5tp_seurat, label = T)
ggsave("./plot_preQC/dimplot_d3to5tp.png",
       width = 8,
       height = 8)
markers <-
  c(
    "Krt5",
    "Krt14",
    "mClover",
    "Acta2",
    "Trp63",
    "Krt18",
    "tdTomato",
    "Kit",
    "Elf5",
    "Ehf",
    "Esr1",
    "Pgr",
    "Foxa1",
    "Sox11",
    "Mki67"
  )
featureplot_d3to5tp_seurat_newmarkers=FeaturePlot(
  d3to5tp_seurat,
  features = markers,
  ncol = 5,
  pt.size = 0.1,
  order = T,
  label = T
)
# dir.create("plot_preQC")
ggsave("./plot_preQC/featureplot_d3to5tp_seurat_newmarkers.png",
       width = 25,
       height = 15)

# pdf(
#   "./plot_preQC/d3to5tp_QCplots.pdf",
#   width = 8.5,
#   height = 7.5,
#   useDingbats = F
# )
VlnPlot(d3to5tp_seurat, features = "nCount_RNA", y.max = 20000) +
  geom_hline(yintercept = 1000) + geom_text(label = "1000", x = 17.25, y =
                                              1000) +
  geom_hline(yintercept = 1250) + geom_text(label = "1250", x = 17.25, y =
                                              1250) +
  geom_hline(yintercept = 1500) + geom_text(label = "1500", x = 17.25, y =
                                              1500)#19 beyound y.max
DimPlot(
  d3to5tp_seurat,
  cells.highlight = WhichCells(d3to5tp_seurat, expression = nCount_RNA <= 1500),
  pt.size = 0.1,
  raster = T
) + ggtitle("nUMI <= 1500 (capped at 20k)")
DimPlot(
  d3to5tp,
  cells.highlight = WhichCells(d3to5tp_seurat, expression = nCount_RNA <= 1500 &
                                 nCount_RNA > 1000),
  pt.size = 0.1,
  raster = T
) + ggtitle("1000 < nUMI <= 1500 (capped at 20k)")
DimPlot(
  d3to5tp_seurat,
  cells.highlight = WhichCells(d3to5tp_seurat, expression = nCount_RNA <= 1500 &
                                 nCount_RNA > 1000 &
                                 nFeature_RNA > 750),
  pt.size = 0.1,
  raster = T
) + ggtitle("1000 < nUMI <= 1500 & nGene > 750(capped at 20k)")# just a few. so nUMI >1k should be fine
## cluter 2 and 8 are low nUMI==> dead cells
## nUMI ~~1250 looks good
## 6 and 11 are high, proliferating
## 9 also contains high nUMI cells==> doublets?
VlnPlot(d3to5tp_seurat, features = "nFeature_RNA") + geom_hline(yintercept = 500) +
  geom_text(label = "500", x = 17.25, y = 500) +
  geom_hline(yintercept = 750, linetype = "dashed") + geom_text(label =
                                                                  "750", x = 17.25, y = 750) +
  geom_hline(yintercept = 1000,
             linetype = "dashed",
             colour = "red") + geom_text(label = "1000", x = 17.25, y = 1000)
DimPlot(
  d3to5tp_seurat,
  cells.highlight = WhichCells(d3to5tp_seurat, expression = nFeature_RNA <= 750),
  pt.size = 0.1,
  raster = T
) + ggtitle("nGene <= 750")# similar to nUMI


VlnPlot(d3to5tp_seurat, features = "percent.mt") + geom_hline(yintercept = 5) +
  geom_hline(yintercept = 7.5, linetype = "dashed") + geom_text(label =
                                                                  "7.5%", x = 17.25, y = 7.5) +
  geom_hline(yintercept = 10,
             linetype = "dashed",
             colour = "red") + geom_text(label = "10%", x = 17.25, y = 10) +
  geom_hline(yintercept = 30,
             linetype = "dashed",
             colour = "blue") + geom_text(label = "30%", x = 17.25, y = 30) +
  geom_hline(yintercept = 35,
             linetype = "dashed",
             colour = "blue") + geom_text(label = "35%", x = 17.25, y = 35)

DimPlot(
  d3to5tp_seurat,
  cells.highlight = WhichCells(d3to5tp_seurat, expression = percent.mt > 35),
  pt.size = 0.1,
  raster = T
) + ggtitle("percent.mt reads > 35%")# similar to nUMI

DimPlot(
  d3to5tp_seurat,
  cells.highlight = WhichCells(d3to5tp_seurat, expression = percent.mt > 35 &
                                 nCount_RNA > 1500),
  raster = T,
  pt.size = 0.1
) + ggtitle("percent.mt reads > 35%  & nUMI >1500")# there are nuclei/debris with high # of nUMI

DimPlot(
  d3to5tp_seurat,
  cells.highlight = WhichCells(d3to5tp_seurat, expression = percent.mt > 35 &
                                 nFeature_RNA > 750),
  pt.size = 0.1
) + ggtitle("percent.mt reads > 35% & nGene > 750")# similar to nUMI

DimPlot(
  d3to5tp_seurat,
  cells.highlight = WhichCells(d3to5tp_seurat, expression = percent.mt > 35 &
                                 nFeature_RNA > 750 &
                                 nCount_RNA > 1500),
  pt.size = 0.1
) + ggtitle("percent.mt reads > 35% & nGene > 750 & nCount_RNA > 1500")# similar to nUMI
dev.off()
# it looks like:
## nGene ~~ 750 looks good
## cluster 2 and 8 are mostly less than 500

dir.create("seurat")
save(d3to5tp_seurat, file = "./seurat/d3to5tp_seurat_raw_seurat.RData")

### qc:
## first REMOVE mitocondrial reads and re-normalize and process
genes <- grep("^mt-", rownames(d3to5tp_seurat), invert = T, value = T)
d3to5tp_seurat.qc <- subset(d3to5tp_seurat, features = genes)
d3to5tp_seurat.qc <- processSeurat(d3to5tp_seurat.qc)
head(d3to5tp_seurat.qc$nCount_RNA)
head(d3to5tp_seurat$nCount_RNA)# automatically updated nCounts
DimPlot(d3to5tp_seurat.qc)# larbely unchanged
markers <-
  c(
    "Krt5",
    "Krt14",
    "Trp63",
    "Acta2",
    "Krt8",
    "Elf5",
    "Ehf",
    "Kit",
    "Esr1",
    "Pgr",
    "Foxa1",
    "Mki67"
  )
FeaturePlot(d3to5tp_seurat.qc,
            features = markers,
            ncol = 4,
            label = T)# LOOKS GOOD BASAL TO LUM TRANSITION
FeaturePlot(
  d3to5tp_seurat.qc,
  features = c("percent.mt", "Mki67", "nCount_RNA"),
  ncol = 3,
  label = T
)#CLUSTER 6 IS HIGH with mt-genes

DimPlot(d3to5tp_seurat.qc, label = T)
dir.create("plot_preQC_mtRm")
ggsave(
  "./plot_preQC_mtRm/d3to5tp_seurat.qc_mt-genes-removed_clusters.png",
  width = 8,
  height = 8
)
markers <-
  c(
    "Krt5",
    "Krt14",
    "mClover",
    "Acta2",
    "Trp63",
    "Krt18",
    "tdTomato",
    "Kit",
    "Elf5",
    "Ehf",
    "Esr1",
    "Pgr",
    "Foxa1",
    "Sox11",
    "Mki67"
  )
FeaturePlot(
  d3to5tp_seurat.qc,
  features = markers,
  ncol = 5,
  pt.size = 0.1,
  order = T,
  label = T
)
ggsave(
  "./plot_preQC_mtRm/d3to5tp_seurat.qc_mt-genes-removed_markers.png",
  width = 25,
  height = 15
)

pdf(
  "./plot_preQC_mtRm//d3to5tp_seurat.qc_mt-genes-removed_QCplots.pdf",
  width = 8.5,
  height = 7.5,
  useDingbats = F
)
VlnPlot(d3to5tp_seurat.qc, features = "nCount_RNA", y.max = 20000) +
  geom_hline(yintercept = 1000) + geom_text(label = "1000", x = 12.25, y =
                                              1000) +
  geom_hline(yintercept = 1250) + geom_text(label = "1250", x = 12.25, y =
                                              1250) +
  geom_hline(yintercept = 1500) + geom_text(label = "1500", x = 12.25, y =
                                              1500)# 1k looks ok, 1.5k better, 15k top

DimPlot(
  d3to5tp_seurat.qc,
  cells.highlight = WhichCells(d3to5tp_seurat.qc, expression = nCount_RNA <= 1500),
  pt.size = 0.1,
  raster = T
) + ggtitle("nUMI <= 1500 (capped at 20k)")
DimPlot(
  d3to5tp_seurat.qc,
  cells.highlight = WhichCells(d3to5tp_seurat.qc, expression = nCount_RNA <= 1500 &
                                 nCount_RNA > 1000),
  pt.size = 0.1,
  raster = T
) + ggtitle("1000 < nUMI <= 1500 (capped at 20k)")
DimPlot(
  d3to5tp_seurat.qc,
  cells.highlight = WhichCells(d3to5tp_seurat.qc, expression = nCount_RNA <= 1500 &
                                 nCount_RNA > 1000 &
                                 nFeature_RNA > 750),
  pt.size = 0.1,
  raster = T
) + ggtitle("1000 < nUMI <= 1500 & nGene > 750(capped at 20k)")# just a few. so nUMI >1k should be fine

VlnPlot(d3to5tp_seurat.qc, features = "nFeature_RNA") + geom_hline(yintercept = 500) +
  geom_text(label = "500", x = 12.25, y = 500) +
  geom_hline(yintercept = 750, linetype = "dashed") + geom_text(label =
                                                                  "750", x = 12.25, y = 750) +
  geom_hline(yintercept = 1000,
             linetype = "dashed",
             colour = "red") + geom_text(label = "1000", x = 12.25, y = 1000)#750 looks great, nGene max 6k

DimPlot(
  d3to5tp_seurat.qc,
  cells.highlight = WhichCells(d3to5tp_seurat.qc, expression = nFeature_RNA <= 750),
  pt.size = 0.1,
  raster = T
) + ggtitle("nGene <= 750")# similar to nUMI


DimPlot(
  d3to5tp_seurat.qc,
  cells.highlight = WhichCells(d3to5tp_seurat.qc, expression = nFeature_RNA < 750 |
                                 nCount_RNA < 1000),
  pt.size = 0.1,
  sizes.highlight = 0.1
) + ggtitle("nGene < 750 & nCount_RNA < 1000")# similar to nUMI
DimPlot(
  d3to5tp_seurat.qc,
  cells.highlight = WhichCells(d3to5tp_seurat.qc, expression = nFeature_RNA < 750 |
                                 nCount_RNA < 1250),
  pt.size = 0.1,
  sizes.highlight = 0.1
) + ggtitle("nGene < 750 & nCount_RNA < 1250")# similar to nUMI
dev.off()

VlnPlot(d3to5tp_seurat.qc, features = "nCount_RNA") + geom_hline(yintercept = 15000)#max 15k

d3to5tp_seurat.qc <-
  subset(d3to5tp_seurat.qc, subset = nFeature_RNA >= 750 &
           nCount_RNA > 1000 & nCount_RNA <= 15000)
d3to5tp_seurat.qc#12141 cells
d3to5tp_seurat.qc <- processSeurat(d3to5tp_seurat.qc)
ElbowPlot(d3to5tp_seurat.qc, ndims = 50)#30 is good enough


pct <- d3to5tp_seurat.qc[["pca"]]@stdev / sum(d3to5tp_seurat.qc[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1








DimPlot(d3to5tp_seurat.qc, label = T)
VlnPlot(d3to5tp_seurat.qc, features = "nCount_RNA")

FeaturePlot(
  d3to5tp_seurat.qc,
  features = c("Epcam"),
  label = T,
  repel = T
)# all epi
FeaturePlot(
  d3to5tp_seurat.qc,
  features = c("Ptprc", "Col1a1", "Dcn", "Kdr", "Rgs5", "Ptprc", "Csf1r", "Cd52"),
  label = T,
  repel = T
)
FeaturePlot(
  d3to5tp_seurat.qc,
  features = c("Mmp7"),
  label = T,
  repel = T
)
# cluster 9 is fibroblast
# cluster 10 is potentiall a doublet
DimPlot(d3to5tp_seurat, label = T)
DimPlot(d3to5tp_seurat.qc, cells.highlight = WhichCells(d3to5tp_seurat, idents = 9))# mainly cluster 10 plus some cluter 1
c11 <-
  FindMarkers(
    d3to5tp_seurat.qc,
    ident.1 = 11,
    logfc.threshold = 1,
    only.pos = T
  )
head(c11, 20)# up3k, mesothelial ==> cluster 11
c12 <-
  FindMarkers(
    d3to5tp_seurat.qc,
    ident.1 = 12,
    logfc.threshold = 1,
    only.pos = T
  )
head(c12, 20)# Cdh19, Negr1==>neuronal 12
FeaturePlot(d3to5tp_seurat.qc, features = rownames(c12)[1:6])
GO <-
  enrichr(rownames(c12), databases = "GO_Biological_Process_2021")
head(GO$GO_Biological_Process_2021)[c(1, 2, 4, 9)]
FeaturePlot(d3to5tp_seurat.qc, features = rownames(c12)[1:6])

c13 <-
  FindMarkers(
    d3to5tp_seurat.qc,
    ident.1 = 13,
    logfc.threshold = 1,
    only.pos = T
  )
head(c13, 20)# skin development
FeaturePlot(
  d3to5tp_seurat.qc,
  features = c(
    "Col17a1",
    "Lamb3",
    "Itgb4",
    "Krt14",
    "Lama3",
    "Itga6",
    "Lamc2",
    "Krt5"
  )
)# nonspecific to c13,
FeaturePlot(d3to5tp_seurat.qc, features = rownames(c13)[1:9])
FeaturePlot(d3to5tp_seurat.qc, features = rownames(c13)[21:29])
FeaturePlot(d3to5tp_seurat.qc, features = rownames(c13)[41:49])
GO <-
  enrichr(rownames(c13)[1:30], databases = "GO_Biological_Process_2021")
head(GO$GO_Biological_Process_2021, 10)[, c(1, 2, 4, 9)]# epidermis
FeaturePlot(d3to5tp_seurat.qc, features = c("Mmp13", "Tll1", "Klk7", "Mmp14"))#
FeaturePlot(d3to5tp_seurat.qc, features = c("Mmp13", "Tll1", "Abca12", "Sprr1b"))#

FeaturePlot(d3to5tp_seurat.qc, features = rownames(c13)[1:6])
FeaturePlot(d3to5tp_seurat.qc, features = "mClover", order = T)
FeaturePlot(d3to5tp_seurat.qc, features = "tdTomato", order = T)
VlnPlot(d3to5tp_seurat.qc, features = "Epcam")
# Ndufa4l2, skin marker
# endou, suprabasal keratinocytes
# DSC3, SUPRABASAL KERATINOCYTES, BASAL
# suprabasal keratinocytes

VlnPlot(d3to5tp_seurat.qc, features = c(rownames(c13)[1:6], "Csf1r", "Cd69"))
VlnPlot(gemm, features = c(rownames(c13)[1:6], "Csf1r", "Cd69"))
load("D:/embryonic_skin/Skin_E13_Ge+Fan_Hormany_integration.RData")
VlnPlot(e13.hn, features = c(rownames(c13)[1:12]))# enriched in those suprabasal cells
DimPlot(e13.hn, group.by = "celltype")

# try doublets finder
library(DoubletFinder)
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_d3to5tp_seurat <-
  paramSweep_v3(d3to5tp_seurat.qc, PCs = 1:30, sct = FALSE)
sweep.stats_d3to5tp_seurat <- summarizeSweep(sweep.res.list_d3to5tp_seurat, GT = FALSE)
bcmvn_d3to5tp_seurat <- find.pK(sweep.stats_d3to5tp_seurat)#0.23

# -------------------------------------------------------------------------------------
homotypic.prop <-
  modelHomotypic(d3to5tp_seurat.qc$seurat_clusters)           ## ex: annotations <- d3to5tp_seurat.qc@meta.data$ClusteringResults
nExp_poi <-
  round(0.08 * nrow(d3to5tp_seurat.qc@meta.data))  ## Assuming 8% doublet
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
d3to5tp_seurat.qc <-
  doubletFinder_v3(
    d3to5tp_seurat.qc,
    PCs = 1:30,
    pN = 0.25,
    pK = 0.23,
    nExp = nExp_poi,
    reuse.pANN = FALSE,
    sct = FALSE
  )
d3to5tp_seurat.qc <-
  doubletFinder_v3(
    d3to5tp_seurat.qc,
    PCs = 1:10,
    pN = 0.25,
    pK = 0.23,
    nExp = nExp_poi.adj,
    reuse.pANN = "pANN_0.25_0.23_971",
    sct = FALSE
  )
colnames(d3to5tp_seurat.qc@meta.data)
FeaturePlot(d3to5tp_seurat.qc, features = "pANN_0.25_0.23_971")# not very obvious
DimPlot(d3to5tp_seurat.qc, label = T, group.by = "DF.classifications_0.25_0.23_808")
FeaturePlot(d3to5tp_seurat.qc, features = "Mki67")

DimPlot(d3to5tp_seurat.qc, label = T)
c5 <-
  FindMarkers(
    d3to5tp_seurat.qc,
    ident.1 = 5,
    logfc.threshold = 1,
    only.pos = T
  )
head(c5)# heatshock markers== heat stressed
c10 <-
  FindMarkers(
    d3to5tp_seurat.qc,
    ident.1 = 10,
    logfc.threshold = 0.7,
    only.pos = T
  )
head(c10)# st18 a tumor suppressor, atp1a2
FeaturePlot(d3to5tp_seurat.qc, features = rownames(c10)[1:9])#
FeaturePlot(d3to5tp_seurat.qc, features = c("Pou2f3", "Dclk1", "Alox5", "Siglecf"))#
FeaturePlot(gemm, features = c("Pou2f3", "Dclk1", "Alox5", "Siglecf"))#
## neuron related genes, atp1a2, dscaml1, st18,
## Pou2f3!
FeaturePlot(d3to5tp_seurat.qc, features = rownames(c5)[1:9])# c10 is not heatshocked
FeaturePlot(d3to5tp_seurat.qc, features = c("Atp1a2", "St18", "Adcy5", "C11orf86"))# neuron markers?
FeaturePlot(d3to5tp_seurat.qc, features = c("Trp63", "Sox10", "Sox11", "Kit"))
library(enrichR)
GO.c10 <-
  enrichr(rownames(c10), databases = "GO_Biological_Process_2021")
head(GO.c10$GO_Biological_Process_2021)[c(1, 2, 4, 9)]
#try scrublets
library(reticulate)
use_condaenv("r-reticulate", conda = "/home/mzhibo/anaconda3/envs/r-reticulate")
scr <- import("scrublet")
d3to5tp_seurat.db.scrublet <-
  scr$Scrublet(
    counts_matrix = t(as.matrix(d3to5tp_seurat.qc@assays$RNA@counts)),
    expected_doublet_rate = 0.1,
    random_state = 2021L
  )
result.d3to5tp_seurat.db.scrublet <-
  d3to5tp_seurat.db.scrublet$scrub_doublets(
    min_counts = 2,
    min_cells = 3,
    min_gene_variability_pctl = 85,
    n_prin_comps = 30L
  )
#d3to5tp_seurat.qc$dbscore.scrublet.08 <- d3to5tp_seurat.db.scrublet[[1]]
d3to5tp_seurat.qc$dbscore.scrublet.10 <- d3to5tp_seurat.db.scrublet[[1]]
#d3to5tp_seurat.qc$dbcall.scrublet.08 <- d3to5tp_seurat.db.scrublet$call_doublets(threshold = quantile(result.d3to5tp_seurat.db.scrublet[[1]],1-0.08))
d3to5tp_seurat.qc$dbcall.scrublet.10 <-
  d3to5tp_seurat.db.scrublet$call_doublets(threshold = quantile(result.d3to5tp_seurat.db.scrublet[[1]], 1 -
                                                        0.1))
#estimated very high doublet fraction: 50% this is consistent with the transitioning states

DimPlot(d3to5tp_seurat.qc, group.by = "dbcall.scrublet.08") +
  DimPlot(d3to5tp_seurat.qc, group.by = "dbcall.scrublet.10")
d3to5tp_seurat.qc$dbcall.scrublet.08 <-
  ifelse(d3to5tp_seurat.qc$dbcall.scrublet.8, "Doublets", "Singlet")
d3to5tp_seurat.qc$dbcall.scrublet.08 <-
  factor(d3to5tp_seurat.qc$dbcall.scrublet.8, levels = c("Singlet", "Doublets"))
d3to5tp_seurat.qc$dbcall.scrublet.10 <-
  ifelse(d3to5tp_seurat.qc$dbcall.scrublet.10, "Doublets", "Singlet")
d3to5tp_seurat.qc$dbcall.scrublet.10 <-
  factor(d3to5tp_seurat.qc$dbcall.scrublet.10, levels = c("Singlet", "Doublets"))

VlnPlot(d3to5tp_seurat.qc, features = "nCount_RNA", split.by = "dbcall.scrublet.8")
VlnPlot(d3to5tp_seurat.qc, features = "nFeature_RNA", split.by = "dbcall.scrublet.8")# looks good
save(d3to5tp_seurat.qc, file = "./d3to5tp_seurat_QCed_seurat.RData")

# remove doublets
d3to5tp_seurat.qc2 <- subset(d3to5tp_seurat.qc, subset = dbcall.scrublet.8 == "Singlet")
d3to5tp_seurat.qc2#11221 cells
d3to5tp_seurat.qc2 <- processSeurat(d3to5tp_seurat.qc2)
DimPlot(d3to5tp_seurat.qc2, label = T)
save(d3to5tp_seurat.qc2, file = "./d3to5tp_seurat_QCed_scrublets-db-8pct_rm_seurat.RData")

d3to5tp_seurat.qc3 <-
  subset(d3to5tp_seurat.qc2, idents = c(9, 12, 11), invert = T)# fibroblast, neuronal, mesothelials removed
d3to5tp_seurat.qc3 <- processSeurat(d3to5tp_seurat.qc3)
DimPlot(d3to5tp_seurat.qc3, label = T)
VlnPlot(d3to5tp_seurat.qc3, features = "nCount_RNA")

load("/mnt/f/sourceData/cc.genes.updated2019_convertedFromHuman.RData")
d3to5tp_seurat.qc3 <-
  CellCycleScoring(d3to5tp_seurat.qc3, s.features = s.genes, g2m.features = g2m.genes)
save(d3to5tp_seurat.qc3, file = "./d3to5tp_seurat_QCed_final_epithelium_seurat_10854cells.RData")

# TRY REGRESSOUT CELLCYCLE
d3to5tp_seurat.qc3.cycrg <-
  ScaleData(
    d3to5tp_seurat.qc3,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = rownames(d3to5tp_seurat.qc3)
  )
d3to5tp_seurat.qc3.cycrg <-
  RunPCA(d3to5tp_seurat.qc3.cycrg) %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims =
                                                                      1:30)
d3to5tp_seurat.qc3.cycrg$oldCluster <- d3to5tp_seurat.qc3.cycrg$seurat_clusters
d3to5tp_seurat.qc3.cycrg <- FindClusters(d3to5tp_seurat.qc3.cycrg, resolution = 0.5)
ElbowPlot(d3to5tp_seurat.qc3.cycrg, ndims = 50)#30 is good
DimPlot(d3to5tp_seurat.qc3.cycrg, label = T)
FeaturePlot(d3to5tp_seurat.qc3.cycrg, features = markers)
save(d3to5tp_seurat.qc3.cycrg, file = "./d3to5tp_seurat_QCed_final_epithelium_seurat_10854cells_cellcycle_regressed.RData")

d3to5tp_seurat.qc4 <-
  subset(d3to5tp_seurat.qc2, idents = c(9, 12, 11, 13), invert = T)# fibroblast, neuronal, mesothelials, skin-superabasal removed
library(future)
plan("multiprocess", workers = 12)
options(future.globals.maxSize = 80 * 1024 ^ 3)
d3to5tp_seurat.qc4 <-
  CellCycleScoring(d3to5tp_seurat.qc4, s.features = s.genes, g2m.features = g2m.genes)
d3to5tp_seurat.qc4 <- processSeurat(d3to5tp_seurat.qc4)
DimPlot(d3to5tp_seurat.qc4, label = T)
save(d3to5tp_seurat.qc4, file = "./d3to5tp_seurat_QCed_final_epithelium_seurat_10832cells.RData")

# TRY REGRESSOUT CELLCYCLE
d3to5tp_seurat.qc4.cycrg <-
  ScaleData(
    d3to5tp_seurat.qc4,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = rownames(d3to5tp_seurat.qc4)
  )
d3to5tp_seurat.qc4.cycrg <-
  RunPCA(d3to5tp_seurat.qc4.cycrg) %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims =
                                                                      1:30)
d3to5tp_seurat.qc4.cycrg <- FindClusters(d3to5tp_seurat.qc4.cycrg, resolution = 0.5)
ElbowPlot(d3to5tp_seurat.qc4.cycrg, ndims = 50)#30 is good
d3to5tp_seurat.qc4.cycrg <- FindClusters(d3to5tp_seurat.qc4.cycrg, resolution = 0.8)
d3to5tp_seurat.qc4.cycrg <- FindClusters(d3to5tp_seurat.qc4.cycrg, resolution = 1)
DimPlot(d3to5tp_seurat.qc4.cycrg, label = T)# looks good
VlnPlot(d3to5tp_seurat.qc4.cycrg, features = c("nCount_RNA", "nFeature_RNA"))#3 AND 10 ARE LOW nUMI, nGenes, #8 is higher
d3to5tp_seurat.qc4.cycrg <-
  RunUMAP(
    d3to5tp_seurat.qc4.cycrg,
    reduction.name = "umap3d",
    reduction.key = "umap3d_",
    dims = 1:30
  )
FeaturePlot(d3to5tp_seurat.qc4.cycrg, features = markers, label = T)

d3to5tp_seurat.qc4.cycrg <-
  BuildClusterTree(
    d3to5tp_seurat.qc4.cycrg,
    dims = 1:30,
    reorder = T,
    reorder.numeric = T
  )#RE-ORDER CLUSTERS
DimPlot(d3to5tp_seurat.qc4.cycrg, label = T)
PlotClusterTree(d3to5tp_seurat.qc4.cycrg)

# annotate
markers.qc4 <-
  FindAllMarkers(
    d3to5tp_seurat.qc4.cycrg,
    min.pct = 0.25,
    logfc.threshold = 0.6,
    only.pos = T
  )# 1.5x difference
write.csv(markers.qc4, file = "./d3to5tp_seurat_QCed_final_epithelium_10832cells_cyc-regressed_12clusters-markers_log2fc0.6_minpct0.25.csv", quote = F)
markers.qc4 <-
  read.csv(
    "./d3to5tp_seurat_QCed_final_epithelium_10832cells_cyc-regressed_12clusters-markers_log2fc0.6_minpct0.25.csv",
    stringsAsFactors = F,
    row.names = 1
  )
head(markers.qc4)
markers.qc4 <- split(markers.qc4, f = markers.qc4$cluster)
DimPlot(d3to5tp_seurat.qc4.cycrg, label = T)
library(enrichR)
head(markers.qc4)
GO.list <-
  lapply(
    lapply(
      markers.qc4,
      FUN = function(x) {
        x$gene
      }
    ),
    enrichr,
    databases = c(
      "GO_Biological_Process_2021",
      "BioCarta_2016",
      "GO_Molecular_Function_2021",
      "WikiPathways_2019_Mouse",
      "WikiPathway_2021_Human"
    )
  )
names(GO.list)


# basal-like clusters
head(markers.qc4$`5`$gene, 100)
head(GO.list$`5`$GO_Biological_Process_2021, 10)[c(1, 2, 4, 9)]# ECM, Muscle contractiion, remodeling, CELL ADHESION
head(GO.list$`5`$BioCarta_2016, 10)[c(1, 2, 4, 9)]#
head(GO.list$`5`$GO_Molecular_Function_2021, 10)[c(1, 2, 4, 9)]# filamin,actin binding, muscle contraction, chemorepellent,
head(GO.list$`5`$WikiPathways_2019_Mouse, 10)[c(1, 2, 4, 9)]# focal adhesinon,adhesion-pik-akt(migration and proliferation, neuron)
head(GO.list$`5`$WikiPathway_2021_Human, 10)[c(1, 2, 4, 9)]# similar

head(markers.qc4$`6`$gene, 100)
head(GO.list$`6`$GO_Biological_Process_2021, 10)[c(1, 2, 4, 9)]# ECM, Muscle contractiion, remodeling, CELL ADHESION, epidermis developmet, axon guidance
head(GO.list$`6`$BioCarta_2016, 10)[c(1, 2, 4, 9)]#
head(GO.list$`6`$GO_Molecular_Function_2021, 10)[c(1, 2, 4, 9)]#
head(GO.list$`6`$WikiPathways_2019_Mouse, 10)[c(1, 2, 4, 9)]# focal adhesinon, spinal cord injury??
head(GO.list$`6`$WikiPathway_2021_Human, 10)[c(1, 2, 4, 9)]# similar

## genes incolved in Hmga2, cellular senscence and chromatin regulatin
## Dclk1, neuronal associated, senscence
## slit2, neuronal guidence and neuronal cell mibration

head(GO.list$`7`$GO_Biological_Process_2021, 10)[c(1, 2, 4, 9)]# similar to cluster 6,5

head(markers.qc4$`7`$gene, 100)
head(GO.list$`7`$GO_Biological_Process_2021, 10)[c(1, 2, 4, 9)]# similar to cluster 6,5




d3to5tp_seurat.qc4.cycrg@reductions$umap_ccnrg <- d3to5tp_seurat.qc4@reductions$umap

save(d3to5tp_seurat.qc4.cycrg, file = "./d3to5tp_seurat_QCed_final_epithelium_seurat_10832cells_cellcycle_regressed.RData")

cellnames.d3to5tp_seurat <- colnames(d3to5tp_seurat.qc4)
cellnames.d3to5tp_seurat.loom <-
  paste0("Multiome_6Dtp:", substr(cellnames.d3to5tp_seurat, 1, 16), "x")
d3to5tp_seurat.qc4.l <- RenameCells(d3to5tp_seurat.qc4, new.names = cellnames.d3to5tp_seurat.loom)

#3 save the umap embedding
write.csv(Embeddings(d3to5tp_seurat.qc4, reduction = "umap"), file = "cell_embeddings_umap-CellCycle-notRegressed_d3to5tp_seurat_loomCellNames.csv")

# export count matrix for pySCENIC analysis
counts <- d3to5tp_seurat.qc4@assays$RNA@counts

counts[1:5, 1:5]
sel <- rowSums(counts > 0)
counts <-
  counts[which(sel > 10),]# only keep genes present in at least 10 cells
write.table(
  counts,
  file = "./pySCENIC/d3to5tp_seurat_QCed_final_epithelium_seurat_10832cells_min10Cells.tsv",
  quote = F,
  col.names = T,
  row.names = T,
  sep = "\t"
)

# label clusters
load(
  "E:/BasalTransplant_Multiome/results/seurat/d3to5tp_seurat_QCed_final_epithelium_seurat_10832cells_cellcycle_regressed.RData"
)
DimPlot(d3to5tp_seurat.qc4.cycrg)
DimPlot(d3to5tp_seurat.qc4, cells.highlight = WhichCells(d3to5tp_seurat.qc4.cycrg, idents = 3))
FeaturePlot(d3to5tp_seurat.qc4, features = c("Mki67", "Stmn1"))
DimPlot(d3to5tp_seurat.qc4, group.by = "Phase")# cluster 3 are not special S phase or G2M phase cells

DimPlot(d3to5tp_seurat.qc4.cycrg, label = T)#
PlotClusterTree(d3to5tp_seurat.qc4.cycrg)
FeaturePlot(
  d3to5tp_seurat.qc4.cycrg,
  features = c("Trp63", "Krt14", "Elf5", "Kit", "Foxa1", "Esr1"),
  label = T,
  ncol = 3,
  pt.size = 0.1
)
clusters <- as.character(Idents(d3to5tp_seurat.qc4.cycrg))
d3to5tp_seurat.qc4.cycrg$celltype_major <-
  ifelse(clusters %in% c(4, 5, 6, 7),
         "Basal-like",
         ifelse(
           clusters %in% c(11),
           "Ba2LP",
           ifelse(
             clusters %in% c(9, 10),
             "LP-like",
             ifelse(
               clusters %in% c(1, 2),
               "ML",
               ifelse(
                 clusters %in% c(8),
                 "LP+ML",
                 ifelse(
                   clusters %in% c(12),
                   "Heat.stressed",
                   ifelse(clusters %in% c(3), "Ba-like.Sox11+", "others")
                 )
               )
             )
           )
         ))
DimPlot(d3to5tp_seurat.qc4.cycrg, label = T, group.by = "celltype_major")
save("d6")

d3to5tp_seurat.qc$seurat_clusters
DimPlot(d3to5tp_seurat.qc, label = T)
DimPlot(d3to5tp_seurat.qc, label = T, group.by = "Phase")
FeaturePlot(d3to5tp_seurat.qc, features = "Mki67")
FeaturePlot(d3to5tp_seurat.qc, features = "Stmn1")
# seperate the prolifrating cells to set higher threshold for AMULET
cells.d3to5tp_seurat.prolif <-
  data.frame(barcode = colnames(d3to5tp_seurat.qc)[which(d3to5tp_seurat.qc$seurat_clusters %in% c(4, 8))],
             is__cell_barcode = 1)
cells.d3to5tp_seurat.nonprolif <-
  data.frame(barcode = colnames(d3to5tp_seurat.qc)[which(!d3to5tp_seurat.qc$seurat_clusters %in% c(4, 8))],
             is__cell_barcode = 1)
write.csv(
  cells.d3to5tp_seurat.prolif,
  "../ArchR/AMULET_doublets/d3to5tp_seurat_scRNAQCed_Cycling-cellBarcodes_1356cells.csv",
  quote = F,
  row.names = F
)
write.csv(
  cells.d3to5tp_seurat.nonprolif,
  "../ArchR/AMULET_doublets/d3to5tp_seurat_scRNAQCed_nonCycling-cellBarcodes_10786cells.csv",
  quote = F,
  row.names = F
)
