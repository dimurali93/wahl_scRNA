processPlots <- function(object, filename,path){
  seu <- object
  # arranging a 5 column grid
  FeaturePlot(seu, features = markers, ncol = 5, label=T)
  ftrPlotFile = paste0(path,filename,"_feature.png")
  ggsave(ftrPlotFile, width = 25, height =15)
  message("Done generating feature plot")
  VlnPlot(seu, features = c("nCount_RNA"))+geom_hline(yintercept = 500)
  nCount_RNA_file = paste0(path,filename,"_nCount_RNA.png")
  ggsave(nCount_RNA_file, width = 25, height =15)
  message("Done generating violin count plot")
  VlnPlot(seu, features = c("nFeature_RNA"))+geom_hline(yintercept = 300)
  nFeature_RNA_file = paste0(path,filename,"_nFeature_RNA.png")
  ggsave(nFeature_RNA_file, width = 25, height =15)
  DimPlot(seu, label = T, 
          cells.highlight = WhichCells(seu,expression = nFeature_RNA >= 300 & nCount_RNA > 500, invert = T), 
          sizes.highlight = 0.1)
  ggsave(paste0(path,filename,".png"), width = 5.5,height = 5)
  message("Done generating feature violin plot")
  return(seu)
}
