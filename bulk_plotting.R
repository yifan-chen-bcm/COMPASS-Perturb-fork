# Importing libraries  
library(tidyverse)
library(readxl)
library(writexl)
library(qs2)
library(Seurat)
library(Mixscale)
library(presto)

# Setting the file path to the 10X files
dir_10x <- "/mnt/mass_storage1/mass_storage_projects/compass/"



# Reading in the processed seurat objects
so <- qs_read(file = paste0(dir_10x, "seurat_objects/2_iNeuron_Mixscale.qs2"))
meta <- so@meta.data



# Visualizations for Mixscale scores per gene (visualizing all gRNAs)
genes <- unique(so$gene)
genes <- genes[genes != "non-targeting"]


for (g in genes) {
  # Subsetting a temporary filtered so
  print(paste0("Plotting: ", g))
  so_temp <- subset(so, gene %in% c("non-targeting", g))
  meta_temp <- so_temp@meta.data
  gRNA_counts <- table(meta_temp$NT_gRNA)
  print(gRNA_counts)
  counts_df <- data.frame(NT_gRNA = names(gRNA_counts), count = as.numeric(gRNA_counts))
  
  # Ridgeplot
  ridge_graph <- RidgePlot(
    so_temp,
    features = "mixscale_score_by_gRNA",
    group.by = "NT_gRNA") + 
    geom_text(
      data = counts_df,
      aes(x = 17.5, y = NT_gRNA, label = paste0("n=", count)),
      inherit.aes = FALSE,
      hjust = 0, vjust = -0.5
    ) +
    labs(title = g) + 
    NoLegend() + 
    coord_cartesian(xlim = c(-8, 20)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  # Saving the Ridgeplot
  ggsave(plot = ridge_graph, 
         filename = paste0("iNeuron_figures/gRNA_ridge/", g, "_mixscale_ridge.pdf"),
         width = 3000, height = 3000, units = "px", dpi = 300)
  
  
  
  # Calculating the max values from violin plot
  expr_max <- FetchData(so_temp, vars = g, slot = "data") %>%
    dplyr::mutate(NT_gRNA = so_temp@meta.data$NT_gRNA) %>%
    dplyr::group_by(NT_gRNA) %>%
    dplyr::summarise(max_expr = max(.data[[g]], na.rm = TRUE))
  
  # Merge with counts
  label_df <- dplyr::left_join(counts_df, expr_max, by = "NT_gRNA") %>%
    dplyr::mutate(ypos = max_expr + 0.1)   # shift labels above violins
  
  # Violin plot
  violin_graph <- VlnPlot(
    so_temp, 
    features = g, 
    group.by = "NT_gRNA", 
    alpha = 0) + 
    geom_text(
      data = label_df,
      aes(x = NT_gRNA, y = ypos, label = paste0("n=", count)),
      inherit.aes = FALSE,
      vjust = 0
    ) +
    geom_boxplot(width = 0.5, outlier.shape = NA) + 
    geom_jitter(width = 0.1, alpha = 0.3, size = 0.5)
  
  # Saving the Violin plot
  ggsave(plot = violin_graph, 
         filename = paste0("iNeuron_figures/gRNA_violin/", g, "_mixscale_violin.pdf"),
         width = 3000, height = 2000, units = "px", dpi = 300)
  
  
  
  # Clearing RAM
  gc()
}




















