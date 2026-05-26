library(tidyverse)
library(ggVennDiagram)
library(openxlsx)
library(ggbeeswarm)
library(fgsea)
library(conflicted)
conflicts_prefer(tidyr::unite)
conflicts_prefer(dplyr::filter)
source("~/Documents/GitHub/rna-seq-host-virus/dge/scripts/Functions.R")
source("~/Documents/GitHub/rna-seq-host-virus/dge/scripts/plot_funcs.R")


base_path <- "~/Library/CloudStorage/GoogleDrive-barbarabitarello@gmail.com/My Drive/BitaLab/BitaLab_not_shared/Research/rna-seq-host-virus/data_and_res_gh_repo/"
subs <- "host"
sampSet <- "76samp" #
res_path <- paste0(base_path, "host/76samp/")
(ext <- paste0("_", sampSet, "_", subs))

## Useful links
#https://biostatsquid.com/gene-set-enrichment-analysis/
#https://biostatsquid.com/easy-gene-set-enrichment-analysis-with-fgsea-old-version/
#https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/
#https://www.gungorbudak.com/blog/2016/05/25/computing-significance-of-overlap/
#https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
#http://www.bioinformatics.cc/article/article-content/268/go-enrichment-analysis-for-non-model-organisms/
# https://biostatsquid.com/fgsea-tutorial-gsea/

## Read in objects ----
bg.genes <- rownames(readRDS(paste0(
  base_path, "host/76samp/txi_76samp_host.rds"
))$counts)

conts <- readRDS(paste0(res_path, "deseq.tb.list.contsOfInterest.rds"))[-c(1:4)]
lapply(names(conts), function(x)
  conts[[x]] |>
    mutate(Ranking = ifelse(
      log2FoldChange < 0, -1, ifelse(log2FoldChange > 0, 1, 0)
    ) * -log10(pvalueRaw)) |>
    mutate(DE = padjIHW <= 0.05 & abs(log2FCshrink_ashr) >= 1) |>
    group_by(locus_tag) |>
    dplyr::slice(1) |>
    #dplyr::select(DE, DIR, Ranking, everything()) |>
    dplyr::select(locus_tag, Ranking) |>
    
    arrange(desc(Ranking)) |>
    write_tsv(
      paste0(
        base_path,
        "data/Annotations-host/STRING-DB-INPUT/",
        x,
        "_sorted.tsv"
      ),
      col_names = F
    ))

conts2 <- lapply(names(conts), function(x)
  conts[[x]] |>
    mutate(Ranking = ifelse(
      log2FoldChange < 0, -1, ifelse(log2FoldChange > 0, 1, 0)
    ) * -log10(pvalueRaw)) |>
    mutate(DE = padjIHW <= 0.05 &
             abs(log2FCshrink_ashr) >= 1) |>
    group_by(locus_tag) |>
    dplyr::slice(1) |>
    #dplyr::select(DE, DIR, Ranking, everything()) |>
    dplyr::select(locus_tag, Ranking, everything()) |>
    arrange(desc(Ranking)))
names(conts2) <- unlist(lapply(conts2, function(x)
  x$contrast[1]))
gc()

## Annotation Files ----
#annot <- readRDS(paste0(base_path, "data/Annotations-host/annot_host4.rds"))
#gene_sets<-paste0(base_path, "data/Annotations-host/gprofiler_full_ehuxleyi.ENSG.gmt") |> gmtPathways()

tempfile1 = tempfile()
read_tsv(
  paste0(
    base_path,
    "data/Annotations-host/gprofiler_full_ehuxleyi.ENSG.gmt"
  ),
  col_names = F
) |>
  unite("X0", X1:X2) |>
  ##|>
  # mutate(X0 = gsub(" ", "_", X0)) |>
  write_tsv(tempfile1)
gene_sets <- tempfile1 |> gmtPathways()
gene_sets <- lapply(gene_sets, function(x)
  x[x != "NA"])

# Degs ----


#names(conts2) #for each

runGSEA <- function(x = 'HHQ_inf_t4vDMSO_inf_t4') {
  #x = 'HHQ_inf_t4vDMSO_inf_t4'
  query1 <- conts2[[x]]
  query2 <- query1$Ranking
  names(query2) <- query1$locus_tag
  #query2<-unique(query2)
  #names(query2)<-unique(query1$locus_tag)
  
  #check
  print(identical(query2 , sort(query2, decreasing = T)))
  print(summary(query2))
  plot(query2)
  
  #ggplot(data.frame(gene_symbol = names(query2)[1:50], ranks = query2[1:50]), aes(gene_symbol, ranks)) +
  #  geom_point() +
  #  theme_classic() +
  #  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  tic("fgsea")
  fgseaRes <- fgsea(
    pathways = gene_sets,
    stats = query2,
    scoreType = "std",
    minSize = 10,
    maxSize = 500,
    nproc = 1,
    gseaParam = 0.5
  )
  toc()
  ##  user  system elapsed
  ## 9.361   5.816   5.302
  fgseaRes |> filter(pval <= 0.05) |> nrow() #295 #this changes every time
  fgseaRes |> filter(padj <= 0.05) |> nrow() #157
  tic("collapse")
  collapsedPathways <- collapsePathways(
    fgseaRes[order(pval)][padj < 0.05],
    stats = query2,
    pathways = gene_sets,
    gseaParam = 0.5
  )
  toc()
  return(collapsedPathways)
}
#sanity check. # https://github.com/alserglab/fgsea/issues/13
res <- vector('list', 10)
for (x in 1:10) {
  cat(x, "\n")
  res[[x]] <- runGSEA()
}

sort(table(unlist(lapply(res, function(x)
  x$mainPathways))))
tibble(GO = unlist(lapply(res, function(x)
  x$mainPathways))) |> group_by(GO) |> tally() |> arrange(desc(n)) |> group_by(n) |> tally() |> arrange(desc(nn))

mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
topPathwaysUp <- fgseaRes |> filter(padj <= 0.05) |> filter(ES > 0) |> arrange(padj) |>
  #slice_head(n = 10) |>
  pull(pathway)
topPathwaysUp <- topPathwaysUp[topPathwaysUp %in% mainPathways]
topPathwaysDown <- fgseaRes |>
  filter(padj <= 0.05) |> filter(ES < 0) |>
  arrange(padj) |>
  #slice_head(n = 10) |>
  pull(pathway)
topPathwaysDown <- topPathwaysDown[topPathwaysDown %in% mainPathways]
(topPathways <- c(topPathwaysUp, rev(topPathwaysDown)))
#run the above agian and compare
#ggVennDiagram(list(topPathways, topPathways2))

p <- plotGseaTable(gene_sets[topPathways],
                   stats = query2,
                   fgseaRes = fgseaRes,
                   gseaParam = 0.5)
p
plotGseaTable(gene_sets[topPathwaysUp], stats  = query2, fgseaRes, gseaParam = 0.5)
pdUp <- plotEnrichmentData(pathway = topPathwaysUp[1],
                           stats = query2,
                           gseaParam = 0.5)
with(
  pdUp,
  ggplot(data = curve) +
    geom_line(aes(x = rank, y = ES), color = "green") +
    geom_ribbon(
      data = stats,
      mapping = aes(
        x = rank,
        ymin = 0,
        ymax = stat / maxAbsStat * (spreadES / 4)
      ),
      fill = "grey"
    ) +
    geom_segment(
      data = ticks,
      mapping = aes(
        x = rank,
        y = -spreadES / 16,
        xend = rank,
        yend = spreadES / 16
      ),
      size = 0.2
    ) +
    geom_hline(
      yintercept = posES,
      colour = "red",
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = negES,
      colour = "red",
      linetype = "dashed"
    ) +
    geom_hline(yintercept = 0, colour = "black") +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "grey92")
    ) +
    labs(x = "rank", y = "enrichment score")
)
#uniprotKEGG<-data.table::fread(paste0(base_path, "data/Annotations-host/annot_uniprot_host.tsv"))

#use uniprot API to convert these
### curl --form 'from="UniProtKB_AC-ID"' \
###--form 'to="UniRef100"' \
###--form 'ids="A0A0M4JKJ3,R1DME7,R1FIW4,R1BX15"' \#except include all genes in uniprotKEGG$Uniprot_Accession
#https://rest.uniprot.org/idmapping/run
##saved as idmapping_2025_08_08.tsv

rosetta <- data.table::fread(paste0(base_path, "data/Annotations-host/idmapping_2025_08_08.tsv"))

rosetta <- rosetta |> dplyr::select(Uniprot_accession = From,
                                    ClusterID_UniRef100 = `Cluster ID`)
rosetta <- rosetta |> mutate(ClusterID_UniRef100 = gsub("UniRef100_", "", ClusterID_UniRef100))

left_join(uniprotKEGG, rosetta) |> write_tsv(paste0(base_path, "data/Annotations-host/annot_uniprot_host2.tsv"))

# Ok now that we have everything, replace names
rosetta <- left_join(uniprotKEGG, rosetta)
rosetta |> filter(is.na(Uniprot_accession))
ehuxGO <- paste0(base_path,
                 "data/Annotations-host/result_2025-7-31-19-26-3.gmt") |> gmtPathways()
goa <- data.table::fread(
  "https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/335831.E_huxleyi_CCMP1516.goa"
) |> dplyr::select(V2, V3, V11)
# convert Uniref100 IDs into something else

#test ---------
i = 10
length(ehuxGO[[i]]) #65
s1 <- ehuxGO[[i]][which(ehuxGO[[i]] %in% rosetta$ClusterID_UniRef100)]
s2 <- ehuxGO[[i]][which(!(ehuxGO[[i]] %in% rosetta$ClusterID_UniRef100))]

part1 <- left_join(tibble(ClusterID_UniRef100 = s1), rosetta)
temp <- left_join(tibble(ClusterID_UniRef100 = ehuxGO[[10]]), rosetta)
temp1 <- temp |> filter(!(is.na(Uniprot_accession)))
temp |> filter(is.na(Uniprot_accession)) |> pull(ClusterID_UniRef100) # only three, looking them up manually
temp2 <- na.omit(
  temp |> filter(is.na(Uniprot_accession))  |>
    mutate(ClusterID_UniRef100 = c("A0A0D3KLH0", "A0A0D3IK94", NA)) |>
    dplyr::select(ClusterID_UniRef100)
)

temp1 <- left_join(temp1 |> dplyr::select(ClusterID_UniRef100) |> distinct(),
                   rosetta)
temp2 <- left_join(temp2 |> dplyr::select(ClusterID_UniRef100) |> distinct(),
                   rosetta)
temp3 <- bind_rows(temp1, temp2)
temp3 |> mutate(gene_names = gsub(" ", ";", gene_names))
temp3 <- temp3 |> mutate(gene_names = gsub(" ", ";", gene_names))

#rename ehuxGO



gsea_test <- fgsea(
  pathways = ehuxGO,
  stats = Rankings,
  scoreType = "std",
  minSize = 10,
  maxSize = 500,
  nproc = 4
)



## Make input files for GSEA with ClusterProfiler

# 1. For each contrast, have gene_ID and shrunken lfc
ldata <- lapply(deseqDegs_lfc1, function(x)
  x |>
    dplyr::select(
      gene_symbol = locus_tag,
      pval = pvalueRaw,
      padj = padjIHW,
      log2fc = log2FCshrink_ashr,
      DIR = DIR
    ))
names(ldata) <- l.namesD
tempF <- lapply(l.namesD, function(x)
  tempfile(pattern = x))
names(tempF) <- l.namesD

ldata2 <- vector('list', length(ldata))
names(ldata2) <- l.namesD
for (cont in l.namesD) {
  ldata2[[cont]] <- list(UP = ldata[[cont]] |> filter(DIR == "UP"),
                         DOWN = ldata[[cont]] |> filter(DIR == "DOWN"))
}
for (cont in l.namesD) {
  write("UP>", file = tempF[[cont]])
  write(gene_symbol)
}
# Then get uniprot identifier for the genes




#virusConts<-readRDS(paste0(base_path, "virus/29samp/tb.list.contsOfInterest.rds"))




### ----- CLusterProfiler
library(AnnotationHub)
ah <- AnnotationHub()
org.Ehux.eg.db <- query(ah, c("OrgDb", "Emiliania"))[[1]]
#orgdb2<-query(ah, c("OrgDb", "Emiliania"))[[2]]
AnnotationDbi::keytypes(org.Ehux.eg.db)
AnnotationDbi::columns(org.Ehux.eg.db)
egid <- keys(orgdb, "ENTREZID")
#egid2 <- keys(orgdb2, "ENTREZID")
AnnotationDbi::select(org.Ehux.eg.db, egid, c("SYMBOL", "GENENAME"), "ENTREZID")
AnnotationDbi::select(org.Ehux.eg.db, egid, c("SYMBOL", "GENENAME"), "ENTREZID") |> dplyr::filter(GENENAME != "hypothetical protein") |> nrow()
#AnnotationDbi::select(orgdb2, egid, c("SYMBOL", "GENENAME"), "ENTREZID") |> dplyr::filter(GENENAME != "hypothetical protein") |> nrow()


tempfile1_in <- read_tsv(tempfile1)

#convert to ENTREZ

geneList <- query1 |> ungroup() |> dplyr::select(Ranking, NCBI_Entrez_geneID) |> na.omit() |> ungroup() |> arrange(desc(Ranking))
Rank <- geneList$Ranking
names(Rank) <- geneList$NCBI_Entrez_geneID
geneList <- Rank
remove(Rank)
res <- clusterProfiler::GSEA(geneList = query2, TERM2GENE = tempfile1_in)
#something wrong here
library(enrichplot)
egoBP <- gseGO(
  geneList = geneList,
  OrgDb        = org.Ehux.eg.db,
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = TRUE
)
dotplot(egoBP)

egoMF <- gseGO(
  geneList = geneList,
  OrgDb        = org.Ehux.eg.db,
  ont          = "MF",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = TRUE
)
dotplot(egoMF)
egoCC <- gseGO(
  geneList = geneList,
  OrgDb        = org.Ehux.eg.db,
  ont          = "CC",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = TRUE
)
dotplot(egoCC)
egoALL <- gseGO(
  geneList = geneList,
  OrgDb        = org.Ehux.eg.db,
  ont          = "ALL",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = TRUE
)
dotplot(egoALL, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free")
egoALL_S <- clusterProfiler::simplify(egoALL)
cnetplot(egoALL_S, foldChange = geneList)
upsetplot(egoALL_S)
ridgeplot(egoALL_S)
ridgeplot(egoALL)
#KEGG
kk <- gseKEGG(query2, nPerm = 10000, organism = "ehx")
ridgeplot(kk)
data.frame(kk)
gseaplot(kk,
         geneSetID = 1,
         by = "runningScore",
         title = kk$Description[1])
gseaplot(kk,
         geneSetID = 1,
         by = "preranked",
         title = kk$Description[1])
dotplotGsea(data = egoALL_S, topn = 15)

library(GseaVis)
gseaNb(object = egoALL_S, geneSetID = 1)