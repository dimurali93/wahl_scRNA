plot_integrated_clusters = function (srat) { 
  ## take an integrated Seurat object, plot distributions over orig.ident
 
  
  
  count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$sample.celltype)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)
  
  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  
  colnames(melt_mtx)[2] <- "dataset"
  colourCount = 15
  GetPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  
  p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
  
  ggplot(mtcars) + 
    geom_histogram(aes(factor(hp), fill=factor(hp))) + 
    scale_fill_manual(values = getPalette(colourCount)) +
    theme(legend.position="bottom") +
    guides(fill=guide_legend(nrow=2))
  
  colourCount = 15
  GetPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  
  # http://rstudio-pubs-static.s3.amazonaws.com/5312_98fc1aba2d5740dd849a5ab797cc2c8d.html
  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
    # scale_fill_brewer(palette = "PiYG") +
    scale_fill_manual(values = GetPalette(colourCount)) +
    theme(legend.position="bottom") +
    guides(fill=guide_legend(nrow=2))
    # ylab("Fraction of cells in each dataset") + xlab("Cluster number") + theme(legend.position="top")
  
  p2 + p1 + plot_layout(widths = c(4,1))
  
  
}
