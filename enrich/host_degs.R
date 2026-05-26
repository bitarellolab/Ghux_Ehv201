library(tidyverse)
library(ggVennDiagram)
library(openxlsx)
library(ggbeeswarm)
library(data.table)
source("~/Documents/GitHub/rna-seq-host-virus/dge/scripts/Functions.R")
source("~/Documents/GitHub/rna-seq-host-virus/dge/scripts/plot_funcs.R")
conflicts_prefer(dplyr::filter)

#https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
#http://www.bioinformatics.cc/article/article-content/268/go-enrichment-analysis-for-non-model-organisms/
#https://seq-jchoi-bio.github.io/docs/RNASeq/HEATMAP_MANUAL/

#tx_to_gene |> filter(!(is.na(`GOs`)) | !(is.na(`GO (GO)`)))

base.path <- "~/Library/CloudStorage/GoogleDrive-barbarabitarello@gmail.com/My Drive/rna-seq-host-virus/data_and_res_gh_repo/"
subs <- "host"
sampSet <- "76samp" #
res_path <- paste0(base.path, "host/76samp/")
(ext <- paste0("_", sampSet, "_", subs))
#degs
res<-readRDS(paste0(res_path, "deseq.tb.list.contsOfInterest.rds"))

#remove sheets that don't contain contrasts 
#names(res)
res<-res[-c(1:4)]
l.names<-names(res)
# End of Deseq X edge -----

###### Deseq Upset plots -----
if(FALSE){
## deseq degs (pIHW<=0.05, |lfc_shrunk| >=1)
deseqDegs_lfc1 <- readRDS(paste0(res_path, "deseq_JustDEGs_p0.05lfc1.rds"))
## edger deqs (padj<=0.05, |lfc| >=0.5
edgeRDegs <- readRDS(paste0(res_path, "edgeR_ContsOfInteres2t_p0.05", ext, ".rds"))$glm
edgeRDegs <- edgeRDegs |> mutate(DIR = ifelse(logFC > 0, "UP", ifelse(logFC <
                                                                     0, "DOWN", NA)), .after = "logFC")
l.names <- unique(edgeRDegs$contrast)

tabEDGE <- edgeRDegs |>
  filter(abs(logFC) >= 1)  |>
  group_by(contrast, DIR, .drop = F) |>
  tally() |>
  pivot_wider(names_from = DIR, values_from = n) |>
  separate(
    contrast,
    into = c("exp", "cntl"),
    sep = "v",
    remove = F
  ) |>
  mutate(DOWN = ifelse(is.na(DOWN) == T, 0, DOWN)) |>
  mutate(UP = ifelse(is.na(UP) == T, 0, UP)) |>
  mutate(Total = DOWN + UP) |>
  mutate(contrast = factor(contrast, levels = names(deseqDegs_lfc1))) |>
  ungroup() |>
  dplyr::select(contrast, exp, cntl, UP, DOWN, Total) |>
  arrange(contrast)
tabEDGE |> gt() |> gtsave(paste0(res_path, "tables/edgeR_DEGsPerContOfInteres2tp0.05lfc1.png"))
tabDeseq <- do.call(rbind, deseqDegs_lfc1) |>
  #filter(contrast %in% l.names) |>
  group_by(contrast, DIR) |>
  tally() |>
  pivot_wider(names_from = DIR, values_from = n) |>
  separate(
    contrast,
    into = c("exp", "cntl"),
    sep = "v",
    remove = F
  ) |>
  mutate(Total = DOWN + UP) |>
  mutate(contrast = factor(contrast, levels = names(deseqDegs_lfc1))) |>
  ungroup() |>
  dplyr::select(contrast, exp, cntl, UP, DOWN, Total) |>
  arrange(contrast)

(
  tabBoth <- left_join(tabDeseq, tabEDGE, by = c("contrast", "exp", "cntl")) |> rename(
    Deseq_UP = UP.x,
    EdgeR_UP = UP.y,
    Deseq_DOWN = DOWN.x,
    EdgeR_DOWN = DOWN.y,
    Deseq_Total = Total.x,
    EdgeR_Total = Total.y
  ) |> mutate(
    EdgeR_UP = ifelse(is.na(EdgeR_UP) == T, 0, EdgeR_UP),
    EdgeR_DOWN = ifelse(is.na(EdgeR_DOWN) == T, 0, EdgeR_DOWN),
    EdgeR_Total = EdgeR_DOWN + EdgeR_UP
  ) |>
    mutate(contrast = factor(contrast, levels = names(deseqDegs_lfc1))) |>
    arrange(contrast)
)

tabBoth |> dplyr::select(-exp) |> gt() |>  gtsave(paste0(res_path, "tables/deseqVedgeR_DEGsp0.05lfc1.png"),
                                                  expand = 20)

edgeRDegs_lfc1 <- edgeRDegs |> filter(abs(logFC) >= 1) |> group_by(contrast) |>
  group_split()#not the best way
l.names <- unlist(lapply(edgeRDegs_lfc1, function(x)
  unique(x$contrast)))
names(edgeRDegs_lfc1) <- l.names



subsetCont <- l.namesD[l.namesD %in% l.names]


deseqDegsNames <- lapply(l.namesD, function(x)
  deseqDegs_lfc1[[x]]$locus_tag)
names(deseqDegsNames) <- l.namesD
edgeRDegsNames <- lapply(subsetCont, function(x)
  edgeRDegs_lfc1[[x]]$locus_tag)
names(edgeRDegsNames) <- subsetCont

degsBoth <- vector('list', length(subsetCont))
names(degsBoth) <- subsetCont


for (cont in l.namesD) {
  if (cont %in% subsetCont) {
    png(paste0(res_path, "figs/deseqVEdge_p0.05lfc1_venn", cont, ".png"))
    degsBoth <- list(edgeR = edgeRDegsNames[[cont]], deseq = deseqDegsNames[[cont]])
    print(ggVennDiagram(degsBoth, force_upset = F, set_size = 1))
    dev.off()
  }
  
  
  
}
}


## We are interes2ted in “the effect of HHQ + virus” for each time point,


# how to get this:

# Effect of HHQ vs DMSO) ----

# For each time point:
#Set A|Exp: HHQ_inf|Cntl: DMSO_inf
#Set C|Exp: HHQ_cntl|Cntl: DMSO_cntl
#Set B|shared between A and C
l.t1names <- l.names[grepl("t1", l.names)]
l.t2names <- l.names[grepl("t2", l.names)]
l.t3names <- l.names[grepl("t3", l.names)]
l.t4names <- l.names[grepl("t4", l.names)]
#make a list with genes for each contrast
#to use later for gsea, make these lists ordered
#which(unlist(lapply(lapply(res, function(x) x |> group_by(contrast) |> summarise(bh = sum(padjBH <=0.05, na.rm = T), ihw = sum(padjIHW<=0.05, na.rm = T))), function(y) y$ihw > y$bh)))
#which(unlist(lapply(lapply(res, function(x) x |> group_by(contrast) |> summarise(bh = sum(padjBH <=0.05, na.rm = T), ihw = sum(padjIHW<=0.05, na.rm = T))), function(y) y$ihw == y$bh)))
#which(unlist(lapply(lapply(res, function(x) x |> group_by(contrast) |> summarise(bh = sum(padjBH <=0.05, na.rm = T), ihw = sum(padjIHW<=0.05, na.rm = T))), function(y) y$ihw < y$bh)))


res2 <- lapply(res, function(x) {
  x |> 
    mutate(Ranking = ifelse(log2FoldChange < 0, -1, 
                               ifelse(log2FoldChange >0, 1, 0))*-log10(pvalueRaw), .after = locus_tag) |> 
    mutate(DE = padjIHW <= 0.05 & abs(log2FCshrink_ashr)>=1, .after = Ranking) |>
    group_by(locus_tag) |>
    dplyr::slice(1) |>
    #dplyr::select(DE, DIR, Ranking, everything()) |> 
    dplyr::select(locus_tag:stat) |>
    arrange(desc(Ranking))
})

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
res2t1<-lapply(res2t1, function(x) x |> filter(DE == T) |> pull(locus_tag))
#custom function
MyUpsetPlot(
  x = res2t1,
  mode = "distinct",
  main = "45 min, DISTINCT",
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
  main = "3 hours, DISTINCT",
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
  main = "8 hours, DISTINCT",
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
  main = "24 hours, DISTINCT",
  file = paste0(res_path, "figs/timePt4upsetDeseqDISTINCT.png")
)
MyUpsetPlot(
  x = res2t4,
  mode = "intersect",
  main = "24 hours, INTERSECT",
  file = paste0(res_path, "figs/timeP4upsetDeseqINTERSECT.png")
)

## Tables ------------------------------
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
# 2. next get effect of virus (inf vs cntl (non-inf)) -----
#For each time point:
#Set D| Exp: HHQ_inf|Cntl: HHQ_cntl)
#Set F| Exp: DMSO_inf|Cntl: DMSO_cntl
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
## Finally Make The Tables ----
annot <- readRDS(paste0(base.path, "data/Annotations-host/annot_", subs, "4.rds"))
Names.AnnotAll<-readxl::excel_sheets(paste0(base.path, "data/Annotations-host/annot-host-ext.xlsx"))
Names.AnnotAll<-Names.AnnotAll[-c(1,2,3,4,9)]
AllAnnot<-vector('list', 4)
names(AllAnnot)<-Names.AnnotAll
for(i in 1:4){
  
try(AllAnnot[[i]]<-readxl::read_xlsx(paste0(base.path, "data/Annotations-host/annot-host-ext.xlsx"), sheet=Names.AnnotAll[i]) )
}
wb1 <- createWorkbook()

# loop over export tibbles to add worksheets

options("openxlsx.maxWidth" = 40)

for (i in 1:length(names(Setst1))) {
  cat(i, "\n")
  addWorksheet(wb1, sheetName = names(Setst1)[i])
  writeDataTable(
    wb1,
    sheet = names(Setst1)[i],
    left_join(tibble(locus_tag = Setst1[[i]]), annot, multiple = "first", by = "locus_tag"),
    tableStyle = "TablestyleMedium2"
  )
  setColWidths(
    wb1,
    sheet = names(Setst1)[i],
    cols = 1:length(Setst1[[i]]),
    widths = "auto"
  )
}

saveWorkbook(
  wb1,
  overwrite = TRUE,
  file = paste0(res_path, "tables/deseq_DEGSets_t1", ext, ".xlsx")
)

wb2 <- createWorkbook()
options("openxlsx.maxWidth" = 40)

for (i in 1:length(names(Setst2))) {
  cat(i, "\n")
  addWorksheet(wb2, sheetName = names(Setst2)[i])
  writeDataTable(
    wb2,
    sheet = names(Setst2)[i],
    left_join(tibble(locus_tag = Setst2[[i]]), annot, multiple = "first", by = "locus_tag"),
    tableStyle = "TablestyleMedium2"
  )
  setColWidths(
    wb2,
    sheet = names(Setst2)[i],
    cols = 1:length(Setst2[[i]]),
    widths = "auto"
  )
}

saveWorkbook(
  wb2,
  overwrite = TRUE,
  file = paste0(res_path, "tables/deseq_DEGSets_t2", ext, ".xlsx")
)

wb3 <- createWorkbook()
options("openxlsx.maxWidth" = 40)

for (i in 1:length(names(Setst3))) {
  cat(i, "\n")
  addWorksheet(wb3, sheetName = names(Setst3)[i])
  writeDataTable(
    wb3,
    sheet = names(Setst3)[i],
    left_join(tibble(locus_tag = Setst3[[i]]), annot, multiple = 'all', by = "locus_tag"),
    tableStyle = "TablestyleMedium2"
  )
  setColWidths(
    wb3,
    sheet = names(Setst3)[i],
    cols = 1:length(Setst3[[i]]),
    widths = "auto"
  )
}

saveWorkbook(
  wb3,
  overwrite = TRUE,
  file = paste0(res_path, "tables/deseq_DEGSets_t3", ext, ".xlsx")
)

wb4 <- createWorkbook()
options("openxlsx.maxWidth" = 40)

for (i in 1:length(names(Setst4))) {
  cat(i, "\n")
  addWorksheet(wb4, sheetName = names(Setst4)[i])
  writeDataTable(
    wb4,
    sheet = names(Setst4)[i],
    left_join(tibble(locus_tag = Setst4[[i]]), annot, multiple = 'first', by = "locus_tag"),
    tableStyle = "TablestyleMedium2"
  )
  setColWidths(
    wb4,
    sheet = names(Setst4)[i],
    cols = 1:length(Setst4[[i]]),
    widths = "auto"
  )
}

saveWorkbook(
  wb4,
  overwrite = TRUE,
  file = paste0(res_path, "tables/deseq_DEGSets_t4", ext, ".xlsx")
)

gc()

sessionInfo()
#gostplot(gostresUP, capped = FALSE, interactive = TRUE)

#gprofiler2::publish_gosttable(gostresUP)

#gostresUP$result |> as_tibble() |> dplyr::select(term_id, term_name, p_value, significant, term_size, query_size, intersection_size, intersection, everything()) |> write_tsv(paste0(res_path, "/tables/SetA_t4_GO_UP.tsv"))

#resGO[["t4"]][["HHQ_inf_t4vDMSO_inf_t4"]][['UP']]$result |> ggplot(aes(x = source, y = -log10(p_value), color = source)) + geom_jitter(data = resGO[["t4"]][["HHQ_inf_t4vDMSO_inf_t4"]][['UP']]$result |> filter(p_value > 0.05), color = "lightgray") + geom_jitter(data = resGO[["t4"]][["HHQ_inf_t4vDMSO_inf_t4"]][['UP']]$result |> filter(p_value <= 0.05))+bb_theme()

