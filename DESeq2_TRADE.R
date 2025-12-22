# Importing libraries
library(tidyverse)
library(readxl)
library(writexl)
library(qs2)
library(Seurat)
library(Mixscale)
library(SingleCellExperiment)
library(DESeq2)
library(TRADEtools)
library(qqplotr)


# Setting the file path to the 10X files
dir_10x <- "/mnt/mass_storage1/mass_storage_projects/compass/"


# Importing the seurat object
so_filt <- qs_read(file = paste0(dir_10x, "seurat_objects/3_iNeuron_Mixscale_filt.qs2"))


# Aggregate (Pseudobulk) raw counts
so_filt2 <- subset(so_filt, target_gene %in% c("SETD1B", "Non-targeting"))
meta <- so_filt2@meta.data
meta <- meta[, c("id", "gRNA")]
meta <- meta %>% tidyr::separate(col = "gRNA", into = c("targeted_gene", "guide_replicate"), sep = "-")
meta$id <- NULL
so_filt2 <- AddMetaData(so_filt2, metadata = meta)
meta <- so_filt2@meta.data

counts_matrix <- GetAssayData(so_filt2, assay = "RNA", slot = "counts")
pseudo_bulk_counts <- AggregateExpression(so_filt2, assays = "RNA", group.by = c("targeted_gene", "guide_replicate"), slot = "counts", return.seurat = FALSE)$RNA
pseudo_bulk_counts <- as.matrix(pseudo_bulk_counts)

sample_metadata <- data.frame(condition = colnames(pseudo_bulk_counts))
rownames(sample_metadata) <- colnames(pseudo_bulk_counts)
sample_metadata <- sample_metadata %>% tidyr::separate(col = condition, into = c("targeted_gene", "guide_replicate"), sep = "_", remove = FALSE)
sample_metadata$targeted_gene <- factor(sample_metadata$targeted_gene, levels = c("nontargeting", "SETD1B"))
stopifnot(all(colnames(pseudo_bulk_counts) == rownames(sample_metadata)))

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = pseudo_bulk_counts, colData = sample_metadata, design = ~targeted_gene)
keep_genes <- rowSums(counts(dds) >= 10) >= 2  #keep genes with >=10 reads in at least 2 gRNA samples
dds <- dds[keep_genes, ]
dds <- DESeq(dds)
dds_res <- results(dds, format = "DataFrame") %>% as.data.frame()
dds_res <- dds_res %>% rownames_to_column(var = "gene")

# TRADE
TRADE_output <- TRADE(mode = "univariate",
                      results1 = dds_res,
                      annot_table = NULL,
                      genes_exclude = NULL,
                      n_sample = NULL)
TRADE_dis_sum <- TRADE_output$distribution_summary %>% as.data.frame()
TRADE_DEGs <- TRADE_output$significant_genes_FDR %>% as.data.frame()




# qq plots
set.seed(0)
smp <- data.frame(norm = rnorm(10000))

p <- TRADE_DEGs$significant_gene_results_FDR.padj
n <- length(p)
expected <- -log10((1:n) / (n + 1))
p_log <- -log10(sort(p))


plot(expected, p_log,
     xlab = "Expected -log10(p)",
     ylab = "Observed -log10(p)",
     pch = 16, cex = 0.5)
abline(0, 1, lwd = 2)


ggqq <- ggplot(data = df, mapping = aes(sample = p)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
ggqq

ggqq <- ggplot(data = smp, mapping = aes(sample = norm)) +
  geom_qq_band(bandType = "ks", mapping = aes(fill = "KS"), alpha = 0.5) +
  geom_qq_band(bandType = "ts", mapping = aes(fill = "TS"), alpha = 0.5) +
  geom_qq_band(bandType = "pointwise", mapping = aes(fill = "Normal"), alpha = 0.5) +
  geom_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  scale_fill_discrete("Bandtype")
ggqq






# DESeq2 DEGs
#Idents(so_filt) <- "target_gene"
#degs_temp <- FindMarkers(so_filt, ident.1 = "KMT2A", ident.2 = "Non-targeting", logfc.threshold = 0, test.use = "DESeq2")




gata1_url <- "https://github.com/ajaynadig/TRADEtools/releases/download/0.99.0/GATA1_K562Essential_DESeq2output.Rdata"
gata1_file <- tempfile(fileext = ".Rdata")
download.file(gata1_url, gata1_file, mode = "wb")
load(gata1_file)
GATA1_K562_results <- results_pseudobulk


Idents(so) <- "NT_gRNA"
guides <- unique(so$NT_gRNA)
guides <- guides[guides != "Non-targeting"]
gRNA_counts <- table(so$NT_gRNA)
n_total <- length(guides)
n_start = 1
full_degs_list <- list()
target_gene_deg_list <- list()

for (gRNA in guides) {
  
  # Check for min cell count
  if (gRNA_counts[gRNA] < 3) {
    print(paste0(gRNA, " has fewer than 3 cells - skipping. ", n_start, "/", n_total))
    n_start <- n_start + 1
    next
  }
  
  # Wilcox differential expression
  degs_temp <- FindMarkers(so, ident.1 = gRNA, ident.2 = "Non-targeting", logfc.threshold = 0, test.use = "wilcox")
  degs_temp <- rownames_to_column(degs_temp, var = "gene")
  
  # DEG regulation direction
  degs_temp <- degs_temp %>% dplyr::mutate(
    regulation = dplyr::case_when(
      p_val_adj < 0.05 & avg_log2FC > log2(1.3) ~ "up",   #30% increase for upregulated DEGs
      p_val_adj < 0.05 & avg_log2FC < log2(0.7) ~ "down"   #30% decrease for downregulated DEGs
    )
  )
  table_degs <- table(degs_temp$regulation)
  degs_temp$up_DEGs <- table_degs["up"]
  degs_temp$down_DEGs <- table_degs["down"]
  degs_temp$total_DEGs <- sum(unique(degs_temp$up_DEGs), unique(degs_temp$down_DEGs), na.rm = TRUE)
  
  # gRNA info and tidying columns
  degs_temp$gRNA <- gRNA
  degs_temp$cell_count <- gRNA_counts[gRNA]
  degs_temp <- degs_temp %>% tidyr::separate(col = "gRNA", into = c("target_gene", "guide_replicate"), sep = "-", remove = FALSE)
  degs_temp <- degs_temp[, c("gene", "gRNA", "target_gene", "guide_replicate", "pct.1", "pct.2", "avg_log2FC", "p_val", "p_val_adj", "regulation", "up_DEGs", "down_DEGs", "total_DEGs", "cell_count")]
  
  # Adding to list
  full_degs_list[[gRNA]] <- degs_temp
  
  # Filtering for target gene degs only
  gene <- gsub("-g[0-9]+$", "", gRNA)
  degs_temp_filt <- degs_temp[degs_temp$gene == gene, ]
  target_gene_deg_list[[gRNA]] <- degs_temp_filt
  
  # Progress status report
  print(paste0("Done with: ", gRNA, ". ", n_start, "/", n_total))
  n_start <- n_start + 1
}

































