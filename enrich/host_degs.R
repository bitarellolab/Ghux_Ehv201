library(tidyverse)
library(ggVennDiagram)
#library(openxlsx)
library(data.table)
source("dge/Functions.R")
source("dge/plot_funcs.R")
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)

#https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
#http://www.bioinformatics.cc/article/article-content/268/go-enrichment-analysis-for-non-model-organisms/
#https://seq-jchoi-bio.github.io/docs/RNASeq/HEATMAP_MANUAL/

# Goals

#1. Compare deseq and edgeR DEGs for host
#2. Make files with DEGs from list of contrasts of interest
#3. Make other lists of DEGs:
#Set A|Exp: HHQ_inf|Cntl: DMSO_inf
#Set C|Exp: HHQ_cntl|Cntl: DMSO_cntl
#Set B|shared between A and C
#Set D| Exp: HHQ_inf|Cntl: HHQ_cntl)
#Set F| Exp: DMSO_inf|Cntl: DMSO_cntl
#Set E| shared between D and F
#Set A': Uniquely A (i.e., set A minus Set B)
#Set D': Uniquely D (i.e., set D minus set E)
#Set G: shared between A' and D'
  
# Read in stuff -----
base_path <- path.expand("~/Documents/Github/Ehux_Ehv201/scratch/")
subs <- "host"
sampSet <- "76samp" #
res_path <- paste0(base_path, "host/76samp/")
(ext <- paste0("_", sampSet, "_", subs))
res<-readRDS(paste0(res_path, "deseq_glm_CountsContsOfInterest", ext, ".rds"))

# Deseq x EdgeR -------
l.names<-names(res)
lfc<-2
padj<-0.05

tabDESeq<-read_tsv("scratch/host/76samp/tables/Table1.tsv")
degs_deseq_lfc2<-lapply(res, function(x){
  j<-x$contrast[1]
  message("Contrast:", j)
  x |> dplyr::filter(padjIHW <= padj) |> 
    dplyr::filter(abs(log2FCshrink_ashr)>=lfc)})
conts<-names(res)

# EdgeR 

  ## edger deqs (padj<=0.05, |lfc| >=0.5
  edgeR<-read_tsv(paste0(res_path, "tables/edgeR_GLM_DiffExpr_AllContrastsp0.05.tsv.gz")) |>
    filter(contrast %in% conts) |> na.omit()
  edgeR<-edgeR |> 
    mutate(DIR = ifelse(logFC > 0, "UP", ifelse(logFC < 0, "DOWN", NA)), .after = "logFC")
  edgeR<-edgeR |>
    mutate(contrast = factor(contrast, levels=conts))
  edgeR |>
    filter(abs(logFC)>=lfc) |>
    mutate(contrast = factor(contrast, levels=conts)) |>
    group_by(contrast,.drop=F) |> tally() |> 
    arrange(desc(n))|>
    gt()
  l.names <- unique(as.character(edgeR$contrast))
  
  edgeRDegs<-edgeR |>
    filter(abs(logFC)>=lfc) |>
    mutate(contrast = factor(contrast, levels = conts))
  
  tabEDGE <- edgeRDegs |>
    group_by(contrast, DIR) |>
    tally() |>
    pivot_wider(names_from = DIR, values_from = n) |>
    separate(
      contrast,
      into = c("exp", "cntl"),
      sep = "v",
      remove = F
    ) |>
    mutate(Downregulated = ifelse(is.na(DOWN) == T, 0, DOWN)) |>
    mutate(Upregulated = ifelse(is.na(UP) == T, 0, UP)) |>
    mutate(Total = Downregulated + Upregulated) |>
    select(-c("UP", "DOWN", "exp", "cntl")) |> arrange(contrast)
  
  tabEDGE |> ungroup() |> gt() |> gtsave(paste0(res_path, "tables/edgeR_DEGsPerContOfInterestp0.05lfc2.png"))
  tabEDGE |> ungroup() |> write_tsv(paste0(res_path, "tables/edgeR_DEGsPerContOfInterestp0.05lfc2.tsv"))
 

  tabBoth <- left_join(tabDESeq |> ungroup(), tabEDGE |> ungroup(), by = c("contrast")) |> 
    rename(
    Deseq_UP = Upregulated.x,
    EdgeR_UP = Upregulated.y,
    Deseq_DOWN = Downregulated.x,
    EdgeR_DOWN = Downregulated.y,
    EdgeR_Total = Total) |>
    mutate(Deseq_Total = Deseq_UP+Deseq_DOWN) |> 
    mutate(
      EdgeR_UP = ifelse(is.na(EdgeR_UP) == T, 0, EdgeR_UP),
      EdgeR_DOWN = ifelse(is.na(EdgeR_DOWN) == T, 0, EdgeR_DOWN),
      EdgeR_Total = EdgeR_DOWN + EdgeR_UP
    ) |>
    mutate(contrast = factor(contrast, levels = names(degs_deseq_lfc2))) |>
    select(-"Unchanged") |>
    arrange(contrast)
  
write_tsv(tabBoth, file = paste0(res_path, "tables/deseq2xedgeR-degs", ext, ".tsv"))

tabBoth |> 
  select(-contains("Total")) |> 
  gt() |> 
  gtsave(paste0(res_path, "tables/deseq2xedgeR-degs", ext,".png"),
                                                                         expand = 20)

edgeRDegs_lfc2 <- edgeRDegs |> group_by(contrast) |>
  group_split()#not the best way
l.names <- unlist(lapply(edgeRDegs_lfc2, function(x)
  unique(x$contrast)))
names(edgeRDegs_lfc2) <- l.names

subsetCont <- names(edgeRDegs_lfc2)

l.namesD<-names(degs_deseq_lfc2)
deseqDegsNames <- lapply(l.namesD, function(x)
  degs_deseq_lfc2[[x]]$locus_tag)
names(deseqDegsNames) <- l.namesD
edgeRDegsNames <- lapply(subsetCont, function(x)
  edgeRDegs_lfc2[[x]]$locus_tag)
names(edgeRDegsNames) <- subsetCont

degsBoth <- vector('list', length(subsetCont))
names(degsBoth) <- subsetCont

for (cont in subsetCont) {
    message("Contrast:", cont)
    png(paste0(res_path, "figs/deseqVEdge_p0.05lfc2_venn_", cont, ".png"))
    degsBoth <- list(edgeR = edgeRDegsNames[[cont]], deseq = deseqDegsNames[[cont]])
    print(ggVennDiagram(degsBoth, force_upset = F, set_size = 2))
    dev.off()
  }
  
for (cont in subsetCont) {
  message("Contrast:", cont)
  png(paste0(res_path, "figs/deseqVEdge_p0.05lfc2_upset_", cont, ".png"))
  degsBoth <- list(edgeR = edgeRDegsNames[[cont]], deseq = deseqDegsNames[[cont]])
  print(ggVennDiagram(degsBoth, force_upset = T, set_size = 2))
  dev.off()
} 
  


# 3) Sets of DEGs -------
#Set A (HHQ+EhV)|Exp: HHQ_inf|Cntl: DMSO_inf
#Set C (HHQ)|Exp: HHQ_cntl|Cntl: DMSO_cntl
#Set B|shared between A and C
#make a list with genes for each contrast
#to use later for gsea, make these lists ordered
#which(unlist(lapply(lapply(res, function(x) x |> group_by(contrast) |> summarise(bh = sum(padjBH <=0.05, na.rm = T), ihw = sum(padjIHW<=0.05, na.rm = T))), function(y) y$ihw > y$bh)))
#which(unlist(lapply(lapply(res, function(x) x |> group_by(contrast) |> summarise(bh = sum(padjBH <=0.05, na.rm = T), ihw = sum(padjIHW<=0.05, na.rm = T))), function(y) y$ihw == y$bh)))
#which(unlist(lapply(lapply(res, function(x) x |> group_by(contrast) |> summarise(bh = sum(padjBH <=0.05, na.rm = T), ihw = sum(padjIHW<=0.05, na.rm = T))), function(y) y$ihw < y$bh)))



res2 <- lapply(res, function(x) {
  x |> 
    mutate(Ranking = ifelse(log2FoldChange < 0, -1, 
                            ifelse(log2FoldChange >0, 1, 0))*-log10(pvalueRaw), .after = locus_tag) |> 
    mutate(DE = padjIHW <= padj & abs(log2FCshrink_ashr)>=lfc, .after = Ranking) |>
    group_by(locus_tag) |>
    dplyr::slice(1) |>
    #dplyr::select(DE, DIR, Ranking, everything()) |> 
    arrange(desc(Ranking)) |>
    dplyr::select(locus_tag:stat) 
})



# separate time points
l.names<-l.namesD
remove(l.namesD)
l.t1names <- l.names[grepl("t1", l.names)]
l.t2names <- l.names[grepl("t2", l.names)]
l.t3names <- l.names[grepl("t3", l.names)]
l.t4names <- l.names[grepl("t4", l.names)]



res2t1 <- res2[l.t1names]
res2t2 <- res2[l.t2names]
res2t3 <- res2[l.t3names]
res2t4 <- res2[l.t4names]
#make a list with genes for each contrast

#https://tomsing1.github.io/blog/posts/upset_plots/
#also interes2ting: https://jokergoo.github.io/InteractiveComplexHeatmap/articles/deseq2_app.html

names(res2t1) <- gsub("_t1", "", names(res2t1))
names(res2t2) <- gsub("_t2", "", names(res2t2))
names(res2t3) <- gsub("_t3", "", names(res2t3))
names(res2t4) <- gsub("_t4", "", names(res2t4))
#make labels prettier
temp <- as_tibble(do.call(rbind, lapply(names(res2t1), function(x)
  str_split_1(x, "v")))) |> mutate(V1 = PrettyTrtVir(V1), V2 = PrettyTrtVir(V2)) |>
  mutate(contrast = paste0(V1, " vs. ", V2))
names(res2t1) <- temp$contrast
names(res2t2) <- temp$contrast
names(res2t3) <- temp$contrast
names(res2t4) <- temp$contrast

# See this for distinct vs intersect https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
res2t1<-lapply(res2t1, function(x) x |> 
                 filter(DE == T) |> pull(locus_tag))
#custom function
conflicts_prefer(ComplexHeatmap::row_order)

MyUpsetPlot(
  x = res2t1,
  mode = "distinct",
  main = "45 min",
  file = paste0(res_path, "figs/timePt1upsetDeseqDISTINCT.png")
)
MyUpsetPlot(
   x = res2t1,
   mode = "intersect",
   main = "45 min, INTERSECT",
   file = paste0(res_path, "figs/timePt1upsetDeseqINTERSECT.png")
 )
res2t2<-lapply(res2t2, function(x) x |> filter(DE == T) |> pull(locus_tag))
MyUpsetPlot(
  x = res2t2,
  mode = "distinct",
  main = "3 hours",
  file = paste0(res_path, "figs/timePt2upsetDeseqDISTINCT.png")
)
MyUpsetPlot(
  x = res2t2,
  mode = "intersect",
  main = "3 hours, INTERSECT",
  file = paste0(res_path, "figs/timeP2upsetDeseqINTERSECT.png")
)
res2t3<-lapply(res2t3, function(x) x |> filter(DE == T) |> pull(locus_tag))
MyUpsetPlot(
  x = res2t3,
  mode = "distinct",
  main = "8 hours",
  file = paste0(res_path, "figs/timePt3upsetDeseqDISTINCT.png")
)
 MyUpsetPlot(
   x = res2t3,
   mode = "intersect",
   main = "8 hours, INTERSECT",
   file = paste0(res_path, "figs/timeP3upsetDeseqINTERSECT.png")
 )
res2t4<-lapply(res2t4, function(x) x |> filter(DE == T) |> pull(locus_tag))
MyUpsetPlot(
  x = res2t4,
  mode = "distinct",
  main = "24 hours",
  file = paste0(res_path, "figs/timePt4upsetDeseqDISTINCT.png")
)
 MyUpsetPlot(
   x = res2t4,
   mode = "intersect",
   main = "24 hours, INTERSECT",
   file = paste0(res_path, "figs/timeP4upsetDeseqINTERSECT.png")
 )

## Tables
#now make tables
# For each time point:
#Set A|Exp: HHQ_inf|Cntl: DMSO_inf
#Set C|Exp: HHQ_cntl|Cntl: DMSO_cntl
#Set B|shared between A and C
Setst1 <- vector('list', length = 9)
names(Setst1) <- c("SetA",
                   "SetB",
                   "SetC",
                   "SetD",
                   "SetE",
                   "SetF",
                   "SetA_prime",
                   "SetD_prime",
                   "SetG")
Setst1$SetA <- res2t1$`HHQ (inf) vs. DMSO (inf)`
Setst1$SetC <- res2t1$`HHQ (cntl) vs. DMSO (cntl)`
Setst1$SetB <- base::intersect(Setst1$SetA, Setst1$SetC)

Setst2 <- vector('list', length = 9)
names(Setst2) <- c("SetA",
                   "SetB",
                   "SetC",
                   "SetD",
                   "SetE",
                   "SetF",
                   "SetA_prime",
                   "SetD_prime",
                   "SetG")
Setst2$SetA <- res2t2$`HHQ (inf) vs. DMSO (inf)`
Setst2$SetC <- res2t2$`HHQ (cntl) vs. DMSO (cntl)`
Setst2$SetB <- base::intersect(Setst2$SetA, Setst2$SetC)

Setst3 <- vector('list', length = 9)
names(Setst3) <- c("SetA",
                   "SetB",
                   "SetC",
                   "SetD",
                   "SetE",
                   "SetF",
                   "SetA_prime",
                   "SetD_prime",
                   "SetG")
Setst3$SetA <- res2t3$`HHQ (inf) vs. DMSO (inf)`
Setst3$SetC <- res2t3$`HHQ (cntl) vs. DMSO (cntl)`
Setst3$SetB <- base::intersect(Setst3$SetA, Setst3$SetC)

Setst4 <- vector('list', length = 9)
names(Setst4) <- c("SetA",
                   "SetB",
                   "SetC",
                   "SetD",
                   "SetE",
                   "SetF",
                   "SetA_prime",
                   "SetD_prime",
                   "SetG")
Setst4$SetA <- res2t4$`HHQ (inf) vs. DMSO (inf)`
Setst4$SetC <- res2t4$`HHQ (cntl) vs. DMSO (cntl)`
Setst4$SetB <- base::intersect(Setst4$SetA, Setst4$SetC)

## 2. next get effect of virus (inf vs cntl (non-inf)) 
#For each time point:
#Set D (HHQ)| Exp: HHQ_inf|Cntl: HHQ_cntl)
#Set F (Ehv)| Exp: DMSO_inf|Cntl: DMSO_cntl
#Set E| shared between D and F
Setst1$SetD <- res2t1$`HHQ (inf) vs. HHQ (cntl)`
Setst1$SetF <- res2t1$`DMSO (inf) vs. DMSO (cntl)`
Setst1$SetE <- base::intersect(Setst1$SetD, Setst1$SetF)

Setst2$SetD <- res2t2$`HHQ (inf) vs. HHQ (cntl)`
Setst2$SetF <- res2t2$`DMSO (inf) vs. DMSO (cntl)`
Setst2$SetE <- base::intersect(Setst2$SetD, Setst2$SetF)

Setst3$SetD <- res2t3$`HHQ (inf) vs. HHQ (cntl)`
Setst3$SetF <- res2t3$`DMSO (inf) vs. DMSO (cntl)`
Setst3$SetE <- base::intersect(Setst3$SetD, Setst3$SetF)

Setst4$SetD <- res2t4$`HHQ (inf) vs. HHQ (cntl)`
Setst4$SetF <- res2t4$`DMSO (inf) vs. DMSO (cntl)`
Setst4$SetE <- base::intersect(Setst4$SetD, Setst4$SetF)

# 3. Lastly
#For each time point:
#Set A': Uniquely A (i.e., set A minus Set B)
#Set D': Uniquely D (i.e., set D minus set E)
#Set G: shared between A' and D'
Setst1$SetA_prime <- base::setdiff(Setst1$SetA, Setst1$SetB)
Setst1$SetD_prime <- base::setdiff(Setst1$SetD, Setst1$SetE)
Setst1$SetG <- base::intersect(Setst1$SetA_prime, Setst1$SetD_prime)

Setst2$SetA_prime <- base::setdiff(Setst2$SetA, Setst2$SetB)
Setst2$SetD_prime <- base::setdiff(Setst2$SetD, Setst2$SetE)
Setst2$SetG <- base::intersect(Setst2$SetA_prime, Setst2$SetD_prime)

Setst3$SetA_prime <- base::setdiff(Setst3$SetA, Setst3$SetB)
Setst3$SetD_prime <- base::setdiff(Setst3$SetD, Setst3$SetE)
Setst3$SetG <- base::intersect(Setst3$SetA_prime, Setst3$SetD_prime)

Setst4$SetA_prime <- base::setdiff(Setst4$SetA, Setst4$SetB)
Setst4$SetD_prime <- base::setdiff(Setst4$SetD, Setst4$SetE)
Setst4$SetG <- base::intersect(Setst4$SetA_prime, Setst4$SetD_prime)


## make upset with these

names(Setst4)
all_ts <- vector('list', 4)
all_ts[[1]] <- list(
  HHQEffect = Setst1[c("SetA", "SetC", "SetB")],
  VirEffect = Setst1[c("SetD", "SetF", "SetE")],
  HHQandVirEffect = Setst1[c("SetA_prime", "SetD_prime", "SetG")]
)
all_ts[[2]] <- list(
  HHQEffect = Setst2[c("SetA", "SetC", "SetB")],
  VirEffect = Setst2[c("SetD", "SetF", "SetE")],
  HHQandVirEffect = Setst2[c("SetA_prime", "SetD_prime", "SetG")]
)
all_ts[[3]] <- list(
  HHQEffect = Setst3[c("SetA", "SetC", "SetB")],
  VirEffect = Setst3[c("SetD", "SetF", "SetE")],
  HHQandVirEffect = Setst3[c("SetA_prime", "SetD_prime", "SetG")]
)
all_ts[[4]] <- list(
  HHQEffect = Setst4[c("SetA", "SetC", "SetB")],
  VirEffect = Setst4[c("SetD", "SetF", "SetE")],
  HHQandVirEffect = Setst4[c("SetA_prime", "SetD_prime", "SetG")]
)

## Upset plots for sets A-G ----
MyUpsetPlot(
  x = all_ts[[1]][['HHQEffect']],
  main = "HHQ effect, 45min, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/HHQeffectUpsett1.png")
)
MyUpsetPlot(
  x = all_ts[[2]][['HHQEffect']],
  main = "HHQ effect, 3h, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/HHQeffectUpsett2.png")
)
MyUpsetPlot(
  x = all_ts[[3]][['HHQEffect']],
  main = "HHQ effect, 8h, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/HHQeffectUpsett3.png")
)
MyUpsetPlot(
  x = all_ts[[4]][['HHQEffect']],
  main = "HHQ effect, 24h, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/HHQeffectUpsett4.png")
)
#
MyUpsetPlot(
  x = all_ts[[1]][['VirEffect']],
  main = "Vir effect, 45min, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/VireffectUpsett1.png")
)
MyUpsetPlot(
  x = all_ts[[2]][['VirEffect']],
  main = "Vir effect, 3h, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/VireffectUpsett2.png")
)
MyUpsetPlot(
  x = all_ts[[3]][['VirEffect']],
  main = "Vir effect, 8h, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/VireffectUpsett3.png")
)
MyUpsetPlot(
  x = all_ts[[4]][['VirEffect']],
  main = "Vir effect, 24h, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/VireffectUpsett4.png")
)

MyUpsetPlot(
  x = all_ts[[1]][['HHQandVirEffect']],
  main = "HHQ and Vir effect, 45min, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/HHQandVireffectUpsett1.png")
)
MyUpsetPlot(
  x = all_ts[[2]][['HHQandVirEffect']],
  main = "HHQ and Vir effect, 3h, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/HHQandVireffectUpsett2.png")
)
MyUpsetPlot(
  x = all_ts[[3]][['HHQandVirEffect']],
  main = "HHQ and Vir effect, 8h, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/HHQandVireffectUpsett3.png")
)
MyUpsetPlot(
  x = all_ts[[4]][['HHQandVirEffect']],
  main = "HHQ and Vir effect, 24h, distinct",
  mode = "distinct",
  file = paste0(res_path, "figs/HHQandVireffectUpsett4.png")
)
## Tables ----
#annot <- readRDS(paste0(base.path, "data/Annotations-host/annot_", subs, "4.rds"))
annot<-readRDS(paste0("data/Annotations-", subs, "/annot-host-ext2-2026-06-28.rds"))
#Names.AnnotAll<-readxl::excel_sheets(paste0("data/Annotations-host/annot-host-ext-2026-06-28.tsv.gz"))
#Names.AnnotAll<-Names.AnnotAll[-c(1,2,3,4,9)]
#AllAnnot<-vector('list', 4)
#names(AllAnnot)<-Names.AnnotAll
#for(i in 1:4){
  
#  try(AllAnnot[[i]]<-readxl::read_xlsx(paste0(base.path, "data/Annotations-host/annot-host-ext.xlsx"), sheet=Names.AnnotAll[i]) )
#}

# loop over export tibbles to add worksheets

#options("openxlsx.maxWidth" = 40)
wb <- xlsx::createWorkbook()

# loop over export tibbles to add worksheets
jgc <- function()
{
  require(rJava)
  gc()
  .jcall("java/lang/System", method = "gc")
}    

temp<-names(Setst1)
#Set A|Exp: HHQ_inf|Cntl: DMSO_inf
#Set C|Exp: HHQ_cntl|Cntl: DMSO_cntl
#Set B|shared between A and C
#Set D| Exp: HHQ_inf|Cntl: HHQ_cntl)
#Set F| Exp: DMSO_inf|Cntl: DMSO_cntl
#Set E| shared between D and F
#Set A': Uniquely A (i.e., set A minus Set B)
#Set D': Uniquely D (i.e., set D minus set E)
#Set G: shared between A' and D'
temp<-tibble(set=names(Setst1), 
             name = c("HHQ_infvDMSO_inf", 
                      "HHQ_infvDMSO_inf & HHQ_infvHHQ_cntl",
                      "HHQ_cntlvDMSO_cntl", 
                      "HHQ_infvHHQ_cntl", 
                      "HHQ_infvHHQ_cntl & DMSO_infvDMSO_cntl", 
                      "DMSO_infvDMSO_cntl",
                      "HHQ_infvDMSO_inf (unique)",
                      "HHQ_infvHHQ_cntl (unique)",
                      "HHQ_infvDMSO_inf (unique) & HHQ_infvHHQ_cntl (unique)"))
message("Creating sheet: legend")
sh<-xlsx::createSheet(wb, "Legend")
message("Adding data frame ", "Legend")
h<-as.data.frame(temp)
xlsx::addDataFrame(h, sh,row.names = F)
for (i in 1:nrow(temp)) {
  gc()
  jgc()
  j<-temp[i,]$set
  j1<-temp[i,]$name
  message("Creating sheet ", j1)
  sh<-xlsx::createSheet(wb, j)
  message("Adding data frame ", j)
  h<-as.data.frame(left_join(tibble(locus_tag = Setst1[[j]]), annot, multiple = "first", by = "locus_tag"))
  xlsx::addDataFrame(h, sh,row.names = F)
 
}

xlsx::saveWorkbook(
  wb,
  file = paste0(res_path, "tables/deseq_DEGSets_t1", ext, ".xlsx")
)

wb2 <- xlsx::createWorkbook()

message("Creating sheet: legend")
sh<-xlsx::createSheet(wb2, "Legend")
message("Adding data frame ", "Legend")
h<-as.data.frame(temp)
xlsx::addDataFrame(h, sh,row.names = F)

for (i in 1:nrow(temp)) {
  gc()
  jgc()
  j<-temp[i,]$set
  j1<-temp[i,]$name
  message("Creating sheet ", j1)
  sh<-xlsx::createSheet(wb2, j)
  message("Adding data frame ", j)
  h<-as.data.frame(left_join(tibble(locus_tag = Setst2[[j]]), annot, multiple = "first", by = "locus_tag"))
  xlsx::addDataFrame(h, sh,row.names = F)
  
}
xlsx::saveWorkbook(wb2,
  file = paste0(res_path, "tables/deseq_DEGSets_t2", ext, ".xlsx")
)

wb3 <- xlsx::createWorkbook()

message("Creating sheet: legend")
sh<-xlsx::createSheet(wb3, "Legend")
message("Adding data frame ", "Legend")
h<-as.data.frame(temp)
xlsx::addDataFrame(h, sh,row.names = F)

for (i in 1:nrow(temp)) {
  gc()
  jgc()
  j<-temp[i,]$set
  j1<-temp[i,]$name
  message("Creating sheet ", j1)
  sh<-xlsx::createSheet(wb3, j)
  message("Adding data frame ", j)
  h<-as.data.frame(left_join(tibble(locus_tag = Setst3[[j]]), annot, multiple = "first", by = "locus_tag"))
  xlsx::addDataFrame(h, sh,row.names = F)
  
}
xlsx::saveWorkbook(wb3,
                   file = paste0(res_path, "tables/deseq_DEGSets_t3", ext, ".xlsx")
)

wb4 <- xlsx::createWorkbook()

message("Creating sheet: legend")
sh<-xlsx::createSheet(wb4, "Legend")
message("Adding data frame ", "Legend")
h<-as.data.frame(temp)
xlsx::addDataFrame(h, sh,row.names = F)

for (i in 1:nrow(temp)) {
  gc()
  jgc()
  j<-temp[i,]$set
  j1<-temp[i,]$name
  message("Creating sheet ", j1)
  sh<-xlsx::createSheet(wb4, j)
  message("Adding data frame ", j)
  h<-as.data.frame(left_join(tibble(locus_tag = Setst4[[j]]), annot, multiple = "first", by = "locus_tag"))
  xlsx::addDataFrame(h, sh,row.names = F)
  
}
xlsx::saveWorkbook(wb4,
                   file = paste0(res_path, "tables/deseq_DEGSets_t4", ext, ".xlsx")
)
gc()


#gostplot(gostresUP, capped = FALSE, interactive = TRUE)

#gprofiler2::publish_gosttable(gostresUP)

#gostresUP$result |> as_tibble() |> dplyr::select(term_id, term_name, p_value, significant, term_size, query_size, intersection_size, intersection, everything()) |> write_tsv(paste0(res_path, "/tables/SetA_t4_GO_UP.tsv"))

#resGO[["t4"]][["HHQ_inf_t4vDMSO_inf_t4"]][['UP']]$result |> ggplot(aes(x = source, y = -log10(p_value), color = source)) + geom_jitter(data = resGO[["t4"]][["HHQ_inf_t4vDMSO_inf_t4"]][['UP']]$result |> filter(p_value > 0.05), color = "lightgray") + geom_jitter(data = resGO[["t4"]][["HHQ_inf_t4vDMSO_inf_t4"]][['UP']]$result |> filter(p_value <= 0.05))+bb_theme()

## End ----
sessionInfo() |>
  capture.output() |>
  writeLines(paste0("logs/host-degs-sampSet-", lubridate::today(), ".txt"))