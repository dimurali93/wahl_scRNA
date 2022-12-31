# ## To install the package from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("dorothea")
# 
# ## To install the development version from the Github repo:
# devtools::install_github("saezlab/dorothea")
library("dorothea")
library(dorothea)
library(ggplot2)
library(dplyr)
# reference: https://saezlab.github.io/dorothea/articles/single_cell_vignette.html
library(dplyr)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)
library(stringr)

setwd("/home/divya/dimurali/wahl/seurat/dorothea")


output_foldername = "/home/divya/dimurali/wahl/seurat/dorothea"
dorothea_regulon_mm <- get(data("dorothea_mm", package = "dorothea"))
regulon <- dorothea_regulon_mm %>%
  dplyr::filter(confidence %in% c("A","B"))
pbmc=scRNA.merge
# d3to5_subset=subset(x = pbmc, subset = orig.ident == "scrna.D3to5")
d3to5_subset <- run_viper(scRNA.merge, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 1, 
                                 verbose = FALSE))
DefaultAssay(object = d3to5_subset) <- "dorothea"
d3to5_subset <- ScaleData(d3to5_subset)
d3to5_subset <- RunPCA(d3to5_subset, features = rownames(d3to5_subset), verbose = FALSE)
d3to5_subset <- FindNeighbors(d3to5_subset, dims = 1:10, verbose = FALSE)
d3to5_subset <- FindClusters(d3to5_subset, resolution = 0.5, verbose = FALSE)
d3to5_subset <- RunUMAP(d3to5_subset, dims = 1:10, umap.method = "uwot", metric = "cosine")
d3to5_subset_markers <- FindAllMarkers(d3to5_subset, only.pos = FALSE, min.pct = 0.25, 
                               logfc.threshold = 0.25, verbose = FALSE)
DimPlot(d3to5_subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(d3to5_subset, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame() %>%
  t()

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(d3to5_subset)), 
                            cell_type = as.character(d3to5_subset@meta.data$label_scRNA),
                            stringsAsFactors = FALSE)
CellsClusters$cell <- str_replace_all(CellsClusters$cell, '#', '.')
CellsClusters$cell <- str_replace_all(CellsClusters$cell, '-', '.')

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))


# For visualization purposes, we select the 20 most variable TFs across clusters according to our scores.

## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(200, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
palette_length = 100

my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "DoRothEA All tp merged", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 


