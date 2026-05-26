#if (!requireNamespace("BiocManager", quietly=TRUE))
#install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("topGO")

library(topGO)
library(ALL)

base_path <-  "~/Downloads/"
res_path <- paste0(base_path, "temp/76samp/")
source("~/Documents/GitHub/rna-seq-host-virus/dge/scripts/Functions.R")
source("~/Documents/GitHub/rna-seq-host-virus/dge/scripts/plot_funcs.R")


##--pg
#https://bioc.r-universe.dev/articles/topGO/topGO_manual.html
library(topGO)
library(ALL)
#BiocManager::install(affyLib)
library("hgu95av2.db")
data(ALL)

#custom annotations

fil<-"/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/topGO/examples/geneid2go.map"
fil2<-"~/Downloads/temp/13794008/ehux_blast2go_KEGG_merged-mod.csv"
geneID2GO <- readMappings(fil)
str(head(geneID2GO))