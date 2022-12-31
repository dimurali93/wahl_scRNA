## Integration of samples ##########
# normalize and identify variable features for each dataset independently
snRNA.integrate = scRNA.list
snRNA.integrate <- lapply(
  X = snRNA.integrate,
  FUN = function(x) {
    x <- NormalizeData(x)
    x <-
      FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  }
)
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = snRNA.integrate)
## Perform integration
anchors <-
  FindIntegrationAnchors(object.list = snRNA.integrate, anchor.features = features)
# this command creates an 'integrated' data assay
snRNA.combined <- IntegrateData(anchorset = anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(snRNA.combined) <- "integrated"
Idents(snRNA.combined) <- "seurat_clusters"
res = 0.8
message("Running PCA")
# snRNA.combined = FindVariableFeatures(snRNA.combined,nfeatures=2000)
snRNA.combined <- ScaleData(snRNA.combined)
snRNA.combined <- RunPCA(snRNA.combined)
pct <-
  snRNA.combined[["pca"]]@stdev / sum(snRNA.combined[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <-
  sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
ndims <- min(co1, co2)
print(ndims)
message("ndims identified")
message("Running UMAP")
snRNA.combined <-
  RunUMAP(snRNA.combined, dims = 1:ndims, verbose = FALSE)
message("Clustering")
snRNA.combined <-
  FindNeighbors(snRNA.combined, dims = 1:ndims, verbose = FALSE)
snRNA.combined <-
  FindClusters(snRNA.combined, verbose = FALSE, res = res)

pIntegratedBySample = DimPlot(snRNA.combined, label = T, group.by = "orig.ident") + plot_annotation(title = 'Combined clusters by samples')
ggsave("IntegratedBySample.png")
dev.off()

pIntegratedByCluster = DimPlot(snRNA.combined, label = T) + plot_annotation(title = 'Combined clusters')
ggsave("Integrated_Clusters.png")
dev.off()
pIntegratedBySample + pIntegratedByCluster
ggsave("Integrated_All.png")
dev.off()

DefaultAssay(snRNA.combined) <- "RNA"
dir.create(path)
path = "/home/divya/dimurali/wahl_scRNA/CombinedIntegrated_postQC/"
snRNA.combined <-
  genPostFilteredFeaturePlots(snRNA.combined, "CombinedIntegrated_postQC", markers, path)

DimPlot(snRNA.combined, split.by = "orig.ident", label = T)
ggsave("Integrated_SplitBySamples.png")
dev.off()
DimPlot(snRNA.combined, label = T)
ggsave("Integrated_SplitBySamples.png")
dev.off()
FeaturePlot(snRNA.combined, features = "Dcn", split.by = "orig.ident")
ggsave("Integrated_SplitBySamples_DCN.png")
dev.off()
FeaturePlot(snRNA.combined, features = "Sox11", split.by = "orig.ident")
ggsave("Integrated_SplitBySamples_Sox11.png")
dev.off()
FeaturePlot(snRNA.combined, features = "Mki67", split.by = "orig.ident")
ggsave("Integrated_SplitBySamples_Mki67.png")
dev.off()

VlnPlot(snRNA.combined, features = "Sox11", split.by = "orig.ident")
DotPlot(snRNA.combined, features = markers) + RotatedAxis()
ggsave("MarkersVsClusters.png")
dev.off()

df <-
  data.frame(
    snRNA.combined@meta.data$integrated_snn_res.0.8,
    snRNA.combined@meta.data$orig.ident
  )
rownames(df) = rownames(snRNA.combined@meta.data)
colnames(df) <- c("Cell_Type", "Condition")
df <- df %>% group_by(Condition, Cell_Type) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(Percent = Nb / C * 100)
percentClustersMarkers=data.frame(df)
write.csv(percentClustersMarkers,"percentClustersMarkers.csv")
Idents(snRNA.combined) <- "seurat_clusters"
snRNA.combined.markers <- FindAllMarkers(snRNA.combined, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# takes times upt 1 min
write.csv(file="/home/divya/dimurali/wahl_scRNA/seurat/CombinedIntegrated_postQC/ClusterIdentitymarkers_Integrated.csv",snRNA.combined.markers)

snRNA.combined_clusters <- Idents(snRNA.combined)
snRNA.combined$label_snRNA <- snRNA.combined$seurat_clusters
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
new.cluster.ids <-
  c(
    "Ba",
    "cycling",
    "Ba",
    "Ba",
    "Ba",
    "LP",
    "Ba",
    "Ba",
    "Ba",
    "ML",
    "cycling",
    "Ba",
    "cycling+ML",
    "cycling",
    "Ba"
  )
snRNA.combined$label_snRNA <-
  plyr::mapvalues(x = snRNA.combined$label_snRNA,
                  from = current.cluster.ids,
                  to = new.cluster.ids)
# DimPlot(snRNA.list[[1]], group.by = c("label_snRNA","passRNAQC"))+ plot_annotation(title = 'Adult annotated')
p1=DimPlot(snRNA.combined, label = T) + plot_annotation(title = 'IntegratedClusters')
p2=DimPlot(snRNA.combined, group.by = "label_snRNA", label = T) + plot_annotation(title = 'IntegratedSubtypes')
ggarrange(p1,p2)
ggsave("ClusterAnnotation.png")
dev.off()
snRNA.combined$sample.celltype = paste0(snRNA.combined$label_snRNA, ".", snRNA.combined$orig.ident)
snRNA.combined$sample.celltype = gsub(
  pattern = "scrna.",
  replacement = "",
  x = snRNA.combined$sample.celltype
)

DimPlot(snRNA.combined, group.by = "sample.celltype") + plot_annotation(title = 'Integrated Subtypes')
ggsave(
  "Integrated_annotated_clusters_sample.celltype.png",
  width = 25,
  height = 15
)
snRNA.combined = FindVariableFeatures(snRNA.combined,
                                      selection.method = "vst",
                                      nfeatures = 2000)
VariableFeaturePlot(snRNA.combined)
FeatureScatter(
  snRNA.combined,
  feature1 = "PC_1",
  feature2 = "nCount_RNA",
  group.by = "orig.ident"
)#R=-0.07
# nothing informative
# DoHeatmap(
#   snRNA.combined,
#   features = VariableFeatures(snRNA.combined),
#   cells = 1:100,
#   size = 4,
#   assay = "integrated",
#   group.by = "sample.celltype",
#   angle = 90
# ) + NoLegend()
# ggsave("ClusterAnnotation.png")
# dev.off()
save(snRNA.combined, file = "scRNA_QCed_annotated_integrated.RData")


source("~/dimurali/wahl_scRNA/scripts/functions/custom_seurat_functions.R")
plot_integrated_clusters(snRNA.combined)
ggsave("FrequencyPlot.png")


