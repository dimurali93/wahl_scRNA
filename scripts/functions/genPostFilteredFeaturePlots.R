genPostFilteredFeaturePlots <- function(object, filename, markers1,path) {
  seu = object
  FeaturePlot(seu, features = markers1, ncol = 5, label=T)# heatshocked, fibroblast
  ftrPlotFile = paste0(path,filename,"_feature.png")
  ggsave(ftrPlotFile, width = 25, height =15)
  return(seu)
}
