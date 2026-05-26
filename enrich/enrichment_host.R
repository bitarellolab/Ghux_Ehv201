#Enrichment plots
library(gprofiler2)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(conflicted)
conflicts_prefer(dplyr::filter)

#pdf("test.pdf")
base_path <-  "~/Library/CloudStorage/GoogleDrive-barbarabitarello@gmail.com/My Drive/BitaLab/BitaLab_not_shared/Research/rna-seq-host-virus/data_and_res_gh_repo/"
res_path <- paste0(base_path, "host/76samp/")
source("~/Documents/GitHub/rna-seq-host-virus/dge/scripts/Functions.R")
source("~/Documents/GitHub/rna-seq-host-virus/dge/scripts/plot_funcs.R")
#resGO <- readRDS(file = paste0(res_path, "ResGOKeggEnrich.rds"))
deseqDegs_lfc1 <- readRDS(paste0(res_path, "deseq_JustDEGs_p0.05lfc1.rds"))
deseqDegs_lfc1<-lapply(deseqDegs_lfc1, function(x) x |> group_by(locus_tag) |> dplyr::slice(1)) 
pal = c(
  `GO:MF` = "#3B99B1",
  `GO:BP` = "#9FC095",
  `GO:CC` = "#E8A419",
  KEGG = "#F5191C",
  REAC = "#3366cc",
  WP = "#0099c6",
  TF = "#5574a6",
  MIRNA = "#22aa99",
  HPA = "#6633cc",
  CORUM = "#66aa00",
  HP = "#990099"
)

options(ggrepel.max.overlaps = Inf)

# Run GO

resGO <- vector('list', 4)
names(resGO) <- paste0("t", 1:4)


bg.genes<-rownames(readRDS(paste0(res_path, "txi_76samp_host.rds"))$counts)
#new (05/26/2026)
tibble(id=bg.genes) |> write_tsv("~/Downloads/temp/76samp/bg.genes.txt", col_names = F)

ptime<-tibble(time = 1:4, time2 = c("45min", "3h", "8h", "24h"))
labels<-tibble(contrast = names(deseqDegs_lfc1)) |> 
  separate(contrast, into = c("exp", "cntl"), sep = "v", remove = F) |> 
  separate(exp, into = c("trt", "vir", "time"),sep = "_") |> 
  mutate(time = parse_number(time)) |> 
  separate(cntl, into = c("trt2", "vir2", "time2"), sep = "_") |> 
  dplyr::select(-time2) |> left_join(ptime)
labels<-labels |> 
  mutate(title = paste0(trt, " (", vir, ") x ", trt2, " (", vir2, ")")) |> 
  rename(subtitle = time2)

gc()

# Now get all results and write tables

library(openxlsx)

resGO2 <- vector('list', 4)
names(resGO2) <- paste0("t", 1:4)
for (i in 1:4) {
  resGO2[[i]] <- deseqDegs_lfc1[grepl(paste0("t", i), names(deseqDegs_lfc1))]
  wb1 <- createWorkbook()
  options("openxlsx.maxWidth" = 40)
  
  for (j in names(resGO[[i]])) {
    #names(resGO[[i]][[j]])<-c("UP", "DOWN")
    cat(j, "\n")
    
    tmp2 <- deseqDegs_lfc1[[j]] |>
      dplyr::select(locus_tag, DIR, contrast) |>
      distinct()
    cat(nrow(tmp2), " degs", "\n")
    print(tmp2 |> group_by(DIR) |> tally())
    #sanity check
    tmp2$contrast[1] == j
    tmp3<-tmp2 |> group_by(DIR) |> dplyr::select(locus_tag) |> group_split(.keep = T)
    names(tmp3)<-unlist(lapply(tmp3, function(x) unique(x$DIR)))
    tmp3<-lapply(tmp3, function(x) x$locus_tag)
    
    #sanity check
    lapply(tmp3, function(x) sum(x %in% bg.genes))
    
    resGO2[[i]][[j]] <- gost(
      query = tmp3,
      organism = "ehuxleyi",
      exclude_iea = F,
      user_threshold = 0.05,
      correction_method = "g_SCS",
      sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG"),
      highlight = F,
      evcodes = T,
      domain_scope = "custom_annotated",
      significant = F,
      custom_bg = bg.genes
    )
    tmp4<-resGO2[[i]][[j]]$result |> group_by(query) |> group_split()
    names(tmp4)<-unlist(lapply(tmp4, function(x) unique(x$query)))
    for(y in names(tmp4)){
      addWorksheet(wb1, sheetName = paste0(j, "_", y))
      writeDataTable(
        wb1,
        sheet = paste0(j, "_", y),
        tmp4[[y]] |> dplyr::select(DIR = query, significant, padj = p_value, everything()),
        tableStyle = "TablestyleMedium2"
      )
      setColWidths(
        wb1,
        sheet = paste0(j, "_", y),
        cols = 1:ncol(tmp4[[y]]),
        widths = "auto"
      )
    }
  }
  saveWorkbook(
    wb1,
    overwrite = TRUE,
    file = paste0(res_path, paste0("tables/ResGOKeggEnrich_", "t", i, ".xlsx"))
  )
}

#now for plots
for (i in 1:4) {
  resGO[[i]] <- deseqDegs_lfc1[grepl(paste0("t", i), names(deseqDegs_lfc1))]
  
  for (j in names(resGO[[i]])) {
    #names(resGO[[i]][[j]])<-c("UP", "DOWN")
    cat(j, "\n")
    
    tmp2 <- deseqDegs_lfc1[[j]] |>
      dplyr::select(locus_tag, DIR, contrast) |>
      distinct()
    cat(nrow(tmp2), " degs", "\n")
    print(tmp2 |> group_by(DIR) |> tally())
    #sanity check
    tmp2$contrast[1] == j
    tmp3<-tmp2 |> group_by(DIR) |> dplyr::select(locus_tag) |> group_split(.keep = T)
    names(tmp3)<-unlist(lapply(tmp3, function(x) unique(x$DIR)))
    tmp3<-lapply(tmp3, function(x) x$locus_tag)
   
    #sanity check
    lapply(tmp3, function(x) sum(x %in% bg.genes))
    
    resGO[[i]][[j]] <- gost(
      query = tmp3,
      organism = "ehuxleyi",
      exclude_iea = F,
      user_threshold = 0.05,
      correction_method = "g_SCS",
      sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG"),
      highlight = T,
      evcodes = F,
      domain_scope = "custom_annotated",
      significant = T,
      custom_bg = bg.genes
    )
    if(is.null(resGO[[i]][[j]])){
      resGO[[i]][[j]] <- gost(
        query = tmp3,
        organism = "ehuxleyi",
        exclude_iea = F,
        user_threshold = 0.05,
        correction_method = "g_SCS",
        sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG"),
        highlight = F,
        evcodes = F,
        domain_scope = "custom_annotated",
        significant = F,
        custom_bg = bg.genes
      )  
    }

  }
}
saveRDS(resGO, paste0(res_path, "ResGOKeggEnrich.rds"))
#resGO<-readRDS(paste0(res_path, "ResGOKeggEnrich.rds"))


# Plots
for (i in 1:4){
  for (j in names(resGO[[i]])) {
    cat(j, "\n")
    labs <- labels |> 
      filter(contrast == j)
    (p<-MyGoPlot(resGO[[i]][[j]],pal = pal, title = labs$title, subtitle = labs$subtitle) + guides(color = "none", alpha = "none", size = "none"))
    fig.name <- paste0(res_path, "figs/GO_", gsub("_t(1|2|3|4)", "", labs$contrast), "_t",labs$time, ".pdf")
    ggsave(fig.name, width = 11,
           height = 8.5,
           dpi = 300, device = cairo_pdf)
  }
}
gc()

# Tables


sessionInfo()
#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
#wb1 <- createWorkbook()
#options("openxlsx.maxWidth" = 40)

#gostplot(gostresUP, capped = FALSE, interactive = TRUE)

#gprofiler2::publish_gosttable(gostresUP)

#gostresUP$result |> as_tibble() |> dplyr::select(term_id, term_name, p_value, significant, term_size, query_size, intersection_size, intersection, everything()) |> write_tsv(paste0(res_path, "/tables/SetA_t4_GO_UP.tsv"))

#resGO[["t4"]][["HHQ_inf_t4vDMSO_inf_t4"]][['UP']]$result |> ggplot(aes(x = source, y = -log10(p_value), color = source)) + geom_jitter(data = resGO[["t4"]][["HHQ_inf_t4vDMSO_inf_t4"]][['UP']]$result |> filter(p_value > 0.05), color = "lightgray") + geom_jitter(data = resGO[["t4"]][["HHQ_inf_t4vDMSO_inf_t4"]][['UP']]$result |> filter(p_value <= 0.05))+bb_theme()

