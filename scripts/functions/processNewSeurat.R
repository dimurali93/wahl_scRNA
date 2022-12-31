processNewSeurat <- function(parent.object,  res = 0.5) {
  Idents(parent.object)<- "seurat_clusters"
  print(res)
  seu = parent.object
  message("Running PCA") 
  seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures=2000) %>% ScaleData() %>% RunPCA()
  pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  ndims <- min(co1, co2)
  print(ndims)
  message("ndims identified")
  message("Running UMAP")
  seu <- RunUMAP(seu, dims = 1:ndims, verbose = FALSE)
  message("Clustering")
  seu <- FindNeighbors(seu, dims = 1:ndims, verbose = FALSE)
  seu <- FindClusters(seu, verbose = FALSE, res=res)
  return(seu)
}
