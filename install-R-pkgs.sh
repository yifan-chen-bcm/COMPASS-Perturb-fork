#!/bin/bash

# ==========================================
### Set up R working environment and packages for Perturb-seq analysis
# For AWS EC2 Ubuntu 24

# Yifan Chen with Gemini
# 12-14-2025
# ==========================================


# ==========================================
# PART 1: Linux System Dependencies (SUDO)
# ==========================================
echo ">>> Updating Apt and installing System Dependencies..."
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    gfortran \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libglpk-dev \
    libhdf5-dev \
    cmake \
    git


# ==========================================
# PART 2: R Package Installation
# ==========================================
# We use a Here-Document to pass R commands directly to Rscript
echo ">>> Starting R Package Installation. This may take 30-60 minutes."

sudo Rscript -e '
# 1. Setup Mirror and Options
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 2) # Use 2 cores for faster compile

# 2. Install Package Managers
if (!require("devtools")) install.packages("devtools")
if (!require("BiocManager")) install.packages("BiocManager")

# 3. Install CRAN Packages
# Note: "matrixStats" and "igraph" are often explicitly needed for Seurat/WGCNA
cran_pkgs <- c(
    "tidyverse", "readxl", "writexl", "qs2", 
    "Seurat", "pheatmap", "rstatix", "plyr", 
    "gtools", "cowplot", "patchwork", "WGCNA", 
    "matrixStats", "igraph"
)

# Function to install missing CRAN packages
inst <- cran_pkgs[!cran_pkgs %in% installed.packages()[,"Package"]]
if(length(inst)) install.packages(inst)

# 4. Install Bioconductor Packages
bioc_pkgs <- c(
    "DropletUtils", "SingleCellExperiment", 
    "gprofiler2", "glmGamPoi", "impute", 
    "preprocessCore", "GO.db", "limma"
)

BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

# 5. Install GitHub Packages
library(devtools)

# Mixscale
if (!requireNamespace("Mixscale", quietly = TRUE)) {
    tryCatch(
        install_github("satijalab/Mixscale", upgrade = "never"),
        error = function(e) install_github("longmanz/Mixscale", upgrade = "never")
    )
}

# Presto (Fast Wilcoxon)
if (!requireNamespace("presto", quietly = TRUE)) {
    install_github("immunogenomics/presto", upgrade = "never")
}


print(">>> ALL INSTALLATIONS COMPLETE <<<")
'

echo ">>> Script Finished."