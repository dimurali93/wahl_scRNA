# https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/tf_sc.html

regulons <- get_dorothea(organism='mouse', levels=c('A', 'B'))
regulons

d3to5_subset=subset(x = snRNA.combined, subset = orig.ident == "D3to5")
mat <- as.matrix(d3to5_subset@assays$RNA@data)

# this takes atleast 20mins
# Run wmean
wmean <- run_wmean(mat=mat, net=regulons, .source='source', .target='target',
                  .mor='mor', times = 100, minsize = 5)
wmean
# visualisation of interested
d3to5_subset[['tfswmean']] <- wmean%>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = d3to5_subset) <- "tfswmean"

# Scale the data
d3to5_subset <- ScaleData(d3to5_subset)
d3to5_subset@assays$tfswmean@data <- d3to5_subset@assays$tfswmean@scale.data

p1 <- DimPlot(d3to5_subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by ="label_snRNA") 
p2 <- (FeaturePlot(d3to5_subset, features = c("Sox9")) & 
         scale_colour_gradient2(low = 'Green', mid = 'white', high = 'red')) +
  ggtitle('Jun activity')
DefaultAssay(object = d3to5_subset) <- "RNA"
p3 <- FeaturePlot(d3to5_subset, features = c("Sox9")) + ggtitle('Sox9 expression')
DefaultAssay(object = d3to5_subset) <- "tfswmean"
p1 | p2 | p3


# exploration
n_tfs <- 100
# Extract activities from object as a long dataframe
df <- t(as.matrix(d3to5_subset@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = d3to5_subset@meta.data$label_snRNA) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_wmean_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(t(top_wmean_mat), border_color = NA, color=my_color, breaks = my_breaks) 
