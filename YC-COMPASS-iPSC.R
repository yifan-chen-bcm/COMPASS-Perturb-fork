############################################################
# Yifan Chen                                               
# yifan.chen@bcm.edu                                       
# Date: 12-18-2025                 
############################################################

library(tidyverse)
library(ggplot2)
library(readxl)
library(data.table)

# Project title: COMPASS iPSC Analysis

# Importing library
library(tidyverse)
library(readxl)
library(writexl)
library(qs2)
library(Seurat)
library(Mixscale)
library(DropletUtils)
library(SingleCellExperiment)
library(gprofiler2)
library(presto)
library(pheatmap)
library(rstatix)
library(plyr)
library(gtools)
library(cowplot)
library(patchwork)
library(progress)
library(utils)

# Setting the file path to the 10X files
dir_10x <- "/Users/yfchen/COMPASS-Perturb-fork/data/COM-iN-none-aggr/"

