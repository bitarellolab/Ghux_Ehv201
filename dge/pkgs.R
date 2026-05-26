#Packages to load

# The following initializes usage of Bioc devel
#BiocManager::install(version='devel', quietly=T)
# we need BioC v.3.22 (devel) because of one of the packages below.
#https://support.bioconductor.org/p/9161859/
# Update 12/26/25: M Love says to just update all packages: https://support.bioconductor.org/p/9161859/
#Packages
## Run once in a while afte major upgrade:
pkgs <- rownames(installed.packages())
BiocManager::install(pkgs, type = "source", checkBuilt = TRUE)
##

## routine, list bioconductor packages
bioconductor_pks<-c(
  "SummarizedExperiment",
  "DESeq2",
  "tximeta",
  "tximport",
  "PCAtools",
  "DEGreport",
  "sjmisc",
  "ashr",
  "fgsea",
  "edgeR",
  "vsn",
  "BiocParallel",
  "Glimma"
)
BiocManager::install(bioconductor_pks)
## if not installing any packages then update existing with
BiocManager::install(ask = F)
#load
lapply(bioconductor_pks, require, character.only = TRUE)
#other packages
other_packages <-c(
  "tidyverse",
  "data.table",
  "lubridate",
  "openxlsx",
  "Seurat",
  "conflicted",
  "WGCNA"
)
library(pacman)
p_path
pacman::p_load(char=other_packages, install = F)

  #register(MulticoreParam(4))
  # conflicted::conflicts_prefer(dplyr::filter)
  # conflicts_prefer(dplyr::desc)
  # conflicts_prefer(base::intersect)
  # conflicts_prefer(dplyr::rename)
  # conflicts_prefer(generics::as.factor)
  # conflicts_prefer(dplyr::count)
  # conflicts_prefer(MatrixGenerics::rowVarDiffs)
  # conflicts_prefer(stats::cor)
  # conflicts_prefer(pheatmap::pheatmap)
  library("genefilter")
  #conflicts_prefer(genefilter::rowVars)
  # For tibble decimal points, no rounding
  old <- options(
    pillar.sigfig = 6,
    pillar.print_max = 5,
    pillar.print_min = 5,
    pillar.advice = FALSE
)

