filterDeadCells <- function(object) {
  seu = object
  seu$passRNAQC <- TRUE
  # seu$mitoRatio <- PercentageFeatureSet(object = seu, pattern = "^MT-")
  seu=PercentageFeatureSet(seu,pattern = "^mt-", col.name = "percent.mt", assay = "RNA")
  seu$log10GenesperUMI=log10(seu$nFeature_RNA)/log10(seu$nCount_RNA) 
  message("Length before filtering")
  before = length(seu$nFeature_RNA)
  print(before)
  seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA > 300 & percent.mt < 10 & log10GenesperUMI >0.8)# threashold previously explored
  message("Length after filtering")
  after <- length(seu$nFeature_RNA)
  print(after)
  deadcellCnt <- before - after
  msgStr <- paste0("Number of dead cells = ", deadcellCnt)
  message(msgStr)
  return(seu)  
}
