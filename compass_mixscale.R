# Importing libraries
library(tidyverse)
library(readxl)
library(writexl)
library(qs2)
library(Seurat)
library(DoubletFinder)
library(Mixscale)
library(DropletUtils)
library(SingleCellExperiment)
library(gprofiler2)

# Setting the file path to the 10X files
dir_10x <- "/mnt/mass_storage1/mass_storage_projects/compass/"



# Reading in the processed seurat objects
print("Importing seurat object")
so <- qs_read(file = paste0(dir_10x, "seurat_objects/1_iNeuron_processed.qs2"))



# perturbation scores per gRNA
print("Calculating perturb signatures")
so <- CalcPerturbSig(
  object = so,
  assay = "RNA",
  slot = "data",
  gd.class = "NT_gRNA",  #gRNA
  nt.cell.class = "non-targeting",
  reduction = "pca",
  ndims = 5,
  num.neighbors = 20,
  new.assay.name = "Perturb_by_gRNA",  #perturb scores by gRNA
  split.by = NULL  #to specify the metadata column if multiple biological context (like cell lines exist)
)
gc()



# MixScale Scores per gRNA
print("Generating Mixscale scores per gRNA")
so <- RunMixscale(
  object = so,
  assay = "Perturb_by_gRNA",
  slot = "scale.data",
  labels = "NT_gRNA",  #per gRNA
  nt.class.name = "non-targeting",
  min.de.genes = 2,  #5 (default), 2 (works), 1 (FAILS)
  logfc.threshold = 0.2,  #0.2
  de.assay = "RNA",
  max.de.genes = 100,
  new.class.name = "mixscale_score_by_gRNA",  #by gRNAs
  fine.mode = F,
  verbose = T,
  split.by = NULL
)
gc()



# Saving the seurat object
print("Saving seurat object with Mixscale scores")
qs_save(object = so, file = paste0(dir_10x, "seurat_objects/2_iNeuron_Mixscale.qs2"))

# Finishing the script
print("Script Completed!")
































