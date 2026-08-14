#!/opt/R/bin/Rscript
#Packages to load (run once)

# The following initializes usage of Bioc devel
#BiocManager::install(version='devel', quietly=T)
# we need BioC v.3.22 (devel) because of one of the packages below.
#https://support.bioconductor.org/p/9161859/
# Update 12/26/25: M Love says to just update all packages: https://support.bioconductor.org/p/9161859/
#Packages
## Run once in a while after major upgrade ----
#pkgs <- rownames(installed.packages())
#BiocManager::install(pkgs, type = "source", checkBuilt = TRUE)
##


# deseq version 1.49.2 ------------
if(FALSE){ #run once
sha <- "f0e22e34364cf1bb87c8c1d7ffc39b81b822decc"
remotes::install_github("thelovelab/DESeq2", ref = "f0e22e34364cf1bb87c8c1d7ffc39b81b822decc")
library(DESeq2)
}
# edgeR version 4.7.3 -------------
if(FALSE){ #run once
install.packages("~/Library/R/Rlibs/edgeR_4.7.3/edgeR/", repos = NULL, type="source")
library(edgeR)
}

# PCA tools 2.21
if(FALSE){#run once
remotes::install_github("kevinblighe/PCAtools", ref = "1c1d8a997fd4a2977e950cf6bff6704f6aadd824")
library(PCAtools)
}

# Glimma 2.19.2
if(FALSE){#run once
remotes::install_github("mritchielab/GlimmaV2", ref = "f6dad6be2a1e29f67a3553d1317b909b336e52b4")
#library(Glimma)
}
## routine, list bioconductor packages ----
bioconductor_pks<-c(
  "SummarizedExperiment",
#  "DESeq2",
  #"tximeta",
  "tximport",
 # "PCAtools",
#  "DEGreport",
  "sjmisc",
  "ashr",
  "fgsea",
#  "edgeR",
  "vsn",
  "BiocParallel",
 # "Glimma",
  "genefilter",
  "ComplexHeatmap",
  "IHW",
  "WGCNA",
  "AnnotationHub"
)
# run once:
if(FALSE){
  BiocManager::install(bioconductor_pks)
}
## if not installing any packages then update existing with ----
#BiocManager::install(ask = F)
#load
lapply(bioconductor_pks, require, character.only = TRUE)
require(DESeq2)
require(edgeR)
require(PCAtools)
require(Glimma)

#other packages
other_packages <-c(
  "tidyverse",
  "data.table",
  "lubridate",
  "openxlsx",
  "xlsx",
  "conflicted",
  "gprofiler2",
  "patchwork",
  "renv"
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
  # conflicts_prefer(genefilter::rowVars)

sessionInfo()
renv::init() #initiate renv
# Now activate renv (run once)
renv::activate()

