# Mixscale Training

# Load packages
options(Seurat.object.assay.version = 'v3') 
library(Seurat)
library(ggridges)
library(ggplot2)
library(Mixscale)
library(stringr)


# Load example 10x data
# load in the count matrix downloaded from GSE132080
ct_mat = ReadMtx(mtx = "mixscale_training/GSE132080_10X_matrix.mtx", 
                 cells = "mixscale_training/GSE132080_10X_barcodes.tsv", 
                 features = "mixscale_training/GSE132080_10X_genes.tsv")


# load the meta_data
meta_data = read.csv("mixscale_training/GSE132080_cell_identities.csv")
rownames(meta_data) = meta_data$cell_barcode

# create a seurat object 
seurat_obj = CreateSeuratObject(counts = ct_mat, meta.data = meta_data)
rm(ct_mat, meta_data)

# retrieve the guide information for each cell
txt = seurat_obj$guide_identity
txt2 = str_extract(txt, "^[^_]+")
txt3 = gsub(pattern = "^[^_]+_", replacement = "", txt)
seurat_obj[['gene']] = txt2
seurat_obj[['gRNA_name']] = txt3
seurat_obj[['cell_type']] = "K562"
rm(txt, txt2, txt3)

# remove ambiguous cells 
seurat_obj = subset(seurat_obj, subset = number_of_cells == 1)  # 19594 cells remain
seurat_obj = subset(seurat_obj, subset = guide_identity != '*')  # 19587 cells remain

### Pre-processing 
seurat_obj = NormalizeData(seurat_obj)
seurat_obj = FindVariableFeatures(seurat_obj)
seurat_obj = ScaleData(seurat_obj)
seurat_obj = RunPCA(seurat_obj)


# Check seurat_obj meta data
meta.data <- seurat_obj@meta.data
unique(seurat_obj@meta.data$gene) # how many genes targeted
unique(seurat_obj@meta.data$gRNA_name) # guide RNA names

### MixScale


# 1.Calculate Perturbation signatures 


seurat_obj <- CalcPerturbSig(
  object = seurat_obj, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="gene", 
  nt.cell.class = "neg", # set neg control in the gene slot
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, # Briefly speaking, for each cell we will search for its 20 nearest neighbors from the non-targeted (NT) cells, and then remove all technical variation so that perturbation-specific effect can be revealed
  new.assay.name = "PRTB", 
  split.by = NULL)  # to specify the metadata column if multiple biological contexts (e.g., cell lines) exist

seurat_obj@assays # see an additional assay as "PBTR"

# 2. Calculate MixScale Score

seurat_obj = RunMixscale(
  object = seurat_obj, 
  assay = "PRTB",  # choose perturb 
  slot = "scale.data", 
  labels = "gene", 
  nt.class.name = "neg", 
  min.de.genes = 5, 
  logfc.threshold = 0.2,
  de.assay = "RNA",
  max.de.genes = 100, 
  new.class.name = "mixscale_score", 
  fine.mode = F, 
  verbose = F, 
  split.by = NULL)

# Mixscale score will be stored as mixscale_score" in the meta data
seurat_obj@meta.data$mixscale_score


# 3. Check Perturbation score distribution
# a. Check the distribution of the scores for the perturbations
RidgePlot(
  seurat_obj,
  features = "mixscale_score",
  group.by = "gene") + NoLegend()

# b. Check if the scores correlate with the expression level of the target gene itself
Mixscale_ScatterPlot(object = seurat_obj, 
                     nt.class.name = "neg", 
                     slct.ident = unique(seurat_obj$gene)[unique(seurat_obj$gene) != "neg"][1:10], 
                     nbin = 10, 
                     facet_wrap = "gene") + NoLegend()


# 4. Weighted differential gene expression

# run score-based weighted DE test for 12 selected perturbations. It will return a list of data frames (one for each perturbation)
de_res = Run_wmvRegDE(object = seurat_obj, assay = "RNA", slot = "counts",
                      labels = "gene", nt.class.name = "neg", 
                      PRTB_list = c("GATA1", "GINS1", "MTOR"),
                      logfc.threshold = 0.2, 
                      split.by = NULL)

  
# have a quick look at the DE results, stored as a list of dfs
head(de_res[[1]])

gata1_df <- de_res[["GATA1"]]
# select the top 20 DE genes from the perturbation
top_res = de_res[["GATA1"]][order(de_res[["GATA1"]]$p_weight)[1:20], ]
# order the DE genes based on its log-fold-change
top_DEG = rownames(top_res[order(top_res$log2FC, decreasing = T), ])

# heatmap for the top DE genes. cells ordered by Mixscale scores. The expression level of the target gene will be plotted in the first row.
Mixscale_DoHeatmap(object = seurat_obj, 
                   labels = "gene", 
                   nt.class.name = "neg", 
                   slct.ident = "GATA1", 
                   mixscale.score.name = "mixscale_score",
                   features = c("GATA1", top_DEG), angle = 0, hjust = 0.5) 
