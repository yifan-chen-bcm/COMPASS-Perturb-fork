### Install all packages needed
# Yifan Chen, by Gemini


# 1. Install package managers if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("remotes", quietly = TRUE))
  install.packages("remotes")

# 2. Define list of CRAN and Bioconductor packages
# Note: BiocManager::install() works for both CRAN and BioC packages
cran_bioc_packages <- c(
  "tidyverse",
  "readxl",
  "writexl",
  "qs2",
  "Seurat",
  "DropletUtils",        # Bioconductor
  "SingleCellExperiment", # Bioconductor
  "gprofiler2",
  "pheatmap",
  "rstatix",
  "plyr",
  "gtools",
  "cowplot",
  "patchwork",
  "WGCNA"
)

# 3. Install CRAN/Bioc packages
BiocManager::install(cran_bioc_packages)

# 4. Install GitHub-specific packages
# These are commonly used single-cell tools that are not always on CRAN/BioC
# or require specific versions.

# DoubletFinder (User commented out, but included here just in case)
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

# Presto (Fast differential expression)
remotes::install_github("immunogenomics/presto")

# hdWGCNA (High-dimensional WGCNA)
remotes::install_github("smorabit/hdWGCNA")

# Mixscale
# Assuming this refers to the tool for analyzing mix-seq/perturb-seq data
remotes::install_github("longyu-team/Mixscale")

print("Installation complete!")