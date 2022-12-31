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
w8_subset=subset(x = scRNA.merge, subset = orig.ident == "scrna.W8")
w8_subset <- run_viper(w8_subset, regulon,
                          options = list(method = "scale", minsize = 4, 
                                         eset.filter = FALSE, cores = 1, 
                                         verbose = FALSE))
DefaultAssay(object = w8_subset) <- "dorothea"
w8_subset <- ScaleData(w8_subset)
w8_subset <- RunPCA(w8_subset, features = rownames(w8_subset), verbose = FALSE)
w8_subset <- FindNeighbors(w8_subset, dims = 1:10, verbose = FALSE)
w8_subset <- FindClusters(w8_subset, resolution = 0.5, verbose = FALSE)
w8_subset <- RunUMAP(w8_subset, dims = 1:10, umap.method = "uwot", metric = "cosine")
w8_subset_markers <- FindAllMarkers(w8_subset, only.pos = FALSE, min.pct = 0.25, 
                                       logfc.threshold = 0.25, verbose = FALSE)
DimPlot(w8_subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(w8_subset, reduction = "umap",group.by = "label_scRNA", label = TRUE, pt.size = 0.5) 

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df=data.frame()
viper_scores_df <- GetAssayData(w8_subset, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame() %>%
  t()

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(w8_subset)), 
                            cell_type = as.character(w8_subset@meta.data$label_scRNA),
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
                       main = "DoRothEA w8", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 

