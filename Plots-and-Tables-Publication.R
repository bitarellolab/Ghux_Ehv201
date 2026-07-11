# Plots for Publication
source("dge/Functions.R")
source("dge/plot_funcs.R")
library(conflicted)
library(DESeq2)
library(tidyverse)
conflicts_prefer(dplyr::filter)
showtext::showtext_opts(dpi=300)
## this font is finicky
#systemfonts::fonts_as_import(family = "Roboto Condensed")
#systemfonts::match_font("Roboto Condensed")
#base_path <- "~/Library/CloudStorage/GoogleDrive-barbarabitarello@gmail.com/My Drive/BitaLab_not_shared/Research/rna-seq-host-virus/data_and_res_gh_repo/"
base_path<-"~/Documents/GitHub/Ehux_Ehv201/scratch/"


# Template ----
if(FALSE){ #not run
plot_path = path.expand(base_path)
filepath = paste0(plot_path, "test2")
# SVG sizes are in inches, not pixels
res = 300
svglite::svglite(paste0(filepath, ".svg"),
                 width = 1400 / res,
                 height = 930 / res)
#plot
#p1<-ggplot(mtcars, aes(mpg, disp, colour = hp)) + geom_point() + geom_smooth()
p1 <- ggplot(mtcars, aes(mpg, disp, colour = hp)) + geom_point() + geom_smooth() + bb_theme()
p1
dev.off()
# Render the svg into a png image with rsvg via magick
img <- magick::image_read_svg(paste0(filepath, ".svg"), width = 1080)
magick::image_write(img, paste0(filepath, ".png"))
#
}

# Table 1 ---------
## Read in stuff 
base_path <- path.expand("~/Documents/Github/Ehux_Ehv201/scratch/")
subs <- "host"
sampSet <- "76samp" #
res_path <- paste0(base_path, "host/76samp/")
(ext <- paste0("_", sampSet, "_", subs))
res<-readRDS(paste0(res_path, "deseq_glm_CountsContsOfInterest", ext, ".rds"))

l.names<-names(res)
lfc<-2
padj<-0.05

#tab deseq
res2_deseq<-do.call(rbind,res)
tabDESeq<-res2_deseq |> 
  #filter(padjIHW <= padj) |> 
  dplyr::filter(padjIHW <= padj) |>
  dplyr::filter(abs(log2FCshrink_ashr)>=lfc)  |> 
  group_by(contrast, DIR) |> 
  tally() |> 
  pivot_wider(names_from = DIR, values_from = n) |>
  separate(contrast, into = c("exp", "cntl"), sep = "v",remove = F) 
if(sum(colnames(tabDESeq)=="UP")==1){
  tabDESeq<-tabDESeq |> mutate(UP = ifelse(is.na(UP), 0, UP))
}else{
  tabDESeq<-tabDESeq |> mutate(UP = 0)
}
if(sum(colnames(tabDESeq)=="DOWN")==1){
  tabDESeq<-tabDESeq |> mutate(DOWN = ifelse(is.na(DOWN), 0, DOWN))
}else{
  tabDESeq<-tabDESeq |> mutate(DOWN = 0)
}
tabDESeq<-tabDESeq |> mutate(Total = DOWN + UP) |> select(-c(exp, cntl)) |> arrange(desc(Total)) |>
  left_join(res2_deseq|> group_by(contrast, DEG=padjIHW <= padj & abs(log2FCshrink_ashr)>=lfc) |> tally(name = "unchanged") |> filter(DEG==F) |> select(contrast, unchanged))

tabDESeq<-tabDESeq|>select(contrast, `Downregulated` = DOWN, `Upregulated` = UP, Unchanged = unchanged)
tabDESeq |> write_tsv("publication/Table1.tsv")


# 12/26: replotting with all time points per KW's request by email.
# SOM Fig S4 ---------
remove(list=ls())
gc()
source("dge/Functions.R")
source("dge/plot_funcs.R")
library(ComplexHeatmap)
library(conflicted)
library(DESeq2)
library(tidyverse)
conflicts_prefer(dplyr::filter)
showtext::showtext_opts(dpi=300)
subs <- "virus"
sampSet <- "39samp"
base_path<-"~/Documents/GitHub/Ehux_Ehv201/scratch/"
res_path <- paste0(base_path, subs, "/", sampSet, "/")
(ext <- paste0("_", sampSet, "_", subs))

# read in data
dds <- readRDS(paste0(res_path, "deseq_dds", ext, ".rds"))
rld_blind <- readRDS(paste0(res_path, "rld_blind", ext, ".rds"))
vsd_blind <- readRDS(paste0(res_path, "vsd_blind", ext, ".rds"))
#rld <- readRDS(paste0(res_path, "rld", ext, ".rds"))
#vsd <- readRDS(paste0(res_path, "vsd", ext, ".rds"))

#colors
trt_vir_cols <- c(paired_cols[2], paired_cols[8])


#if only t2-t4:
if(sampSet=="39samp"){ 
  # Order time points correctly for plots
  
  times <- c("45min", "3h", "8h", "24h")
  
  # if t1-t4
  
  ann_colors = list(
    Treatment = c(DMSO = trt_vir_cols[1], HHQ = trt_vir_cols[2]),
    Time = c(
      `45min`=grays[2],
      `3h` = grays[3],
      `8h` = grays[4],
      `24h` = grays[5]
    )
  )

#mat <- assay(rld)
mat <- assay(rld_blind)
#mat <- assay(vsd_blind)
#rownames(mat)<-trans_cts[match(rownames(mat),trans_cts$locus_tag),]$gene
scaled_mat = t(scale(t(mat)))

coldata <- read_tsv(paste0("data/coldata_", sampSet, ".txt"),
                    show_col_types = FALSE)
df2 <- as.data.frame(colData(dds)[, c("treatment", "timePt")])
colnames(df2)[1] <- 'trt'
df2$timePt <- PrepNamesForPlots(timePt = df2$timePt)$timePt
# arrange by time
df2$timePt2 <- parse_number(gsub("h", "", df2$timePt))
df2<-df2 |> mutate(timePt2=ifelse(timePt2==45, 0.75, timePt2))
df2 <- df2 |> arrange(timePt2) |> dplyr::select(-timePt2)
#
df2 <- df2 |> dplyr::select(Treatment = trt, Time = timePt)
df2$Rep<-SplitSampleName(sample = rownames(df2))$rep

mat<-mat[,rownames(df2)]
scaled_mat<-scaled_mat[,rownames(df2)]
matDist<-dist(mat)
hClustMat<-hclust(matDist, method = "complete")
#plot(hClustMat, labels = F)

#abline(h = 5, col = "brown", lwd = 2) #
#abline(h = 10, col = "brown", lwd = 2) #

#dend = as.dendrogram(hClustMat)
#library(dendextend)
br_gn<-colorRampPalette(c4a(palette= "brewer.br_bg"))(50)
pheatmap::pheatmap(scaled_mat,
                   main = "Test",
                   cluster_rows = T,
                   show_rownames = F,
                   cluster_cols = F,
                   color = br_gn,
                   annotation_col = df2[,-3],
                   gaps_col = c(10,19,29), annotation_colors = ann_colors)

#sanity check
#apply(scaled_mat, MARGIN = 1, mean) %>%                          # calculate the mean per row
#hist(., main = "", xlab = "Z-score values", col = "dodgerblue2")  
#make names prettier
newLabs<-SplitSampleName(sample = colnames(scaled_mat)) |> mutate(NewName = paste0("Rep ", rep))


names(trt_vir_cols)<-c("DMSO", "HHQ")

# change params globally
ht_opt(legend_border = "black")
#ht_opt(RESET = TRUE)




remove(l1); l1<-length(br_gn);
col_fun_brgn<-colorRamp2(breaks = c(-3, 0, 3), colors = c(br_gn[1],br_gn[l1/2], br_gn[l1]))

  #   ht<-ComplexHeatmap::Heatmap(
  #   scaled_mat,
  #   col = col_fun_brgn,
  #   column_title = "Samples",
  #   column_title_side = "bottom",
  #   row_title = "Genes",
  #   show_row_names = F,
  #   show_column_names = F,
  #   column_names_gp = gpar(fontsize = 8),
  #   #column_labels = newLabs$NewName,
  #   row_names_gp = gpar(fontsize = 4),
  #   show_row_dend = F,
  #   cluster_columns = F,
  #   cluster_rows = T,
  #   heatmap_legend_param = list(title = "", at = c(-4,-3,-2, -1, 0, 1,2,3,4),legend_height = unit(6, "cm")),
  #   top_annotation = ha, 
  #   row_km_repeats = 100,
  #   column_km_repeats = 100,
  #   split = 3,gap = unit(3, "mm"),
  #   show_heatmap_legend = T
  #   #column_split = 2
  # )
  # draw(ht)
  # htb <- ComplexHeatmap::Heatmap(
  #   scaled_mat,
  #   col = col_fun_brgn,
  #   column_title = "Samples",
  #   column_title_side = "bottom",
  #   row_title = "Genes",
  #   show_row_names = F,
  #   show_column_names = F,
  #   column_names_gp = gpar(fontsize = 8),
  #   #column_labels = newLabs$NewName,
  #   row_names_gp = gpar(fontsize = 4),
  #   show_row_dend = F,
  #   cluster_columns = T,
  #   cluster_rows = T,
  #   heatmap_legend_param = list(
  #     title = "",
  #     at = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
  #     legend_height = unit(6, "cm")
  #   ),
  #   top_annotation = ha,
  #   row_km_repeats = 100,
  #   column_km_repeats = 100,
  #   split = 3,
  #   gap = unit(3, "mm"),
  #   column_split = 3,
  #   show_heatmap_legend = T
  #   #column_split = 2
  # )
  # 
  ha1<-HeatmapAnnotation(
    Treatment =df2 |> dplyr::filter(Time == "45min") |> pull(Treatment), Time = df2 |> dplyr::filter(Time == "45min") |> pull(Time),
    col = ann_colors,
    annotation_name_side = "left",which = "column",
    show_legend = c("Treatment" = F, "Time" = F))
  
  ht1<-ComplexHeatmap::Heatmap(
    scaled_mat[,colnames(scaled_mat)[grepl("t1", colnames(scaled_mat))]],
    col = col_fun_brgn,
    column_title = "",
    column_title_side = "bottom",
    row_title = "Genes",
    show_row_names = F,
    show_column_names = F,
    column_names_gp = gpar(fontsize = 8),
    #column_labels = newLabs$NewName,
    row_names_gp = gpar(fontsize = 4),
    show_row_dend = T,
    cluster_columns = F,
    cluster_rows = T,
    heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")),
    top_annotation = ha1, 
    row_km_repeats = 100,
    column_km_repeats = 100,
    split = 3,
    gap = unit(3, "mm"),
    show_heatmap_legend = F
    #column_split = 2
  )
  
  draw(ht1)
  
  ha2<-HeatmapAnnotation(
    Treatment = df2 |> 
      dplyr::filter(Time == "3h") |> 
      pull(Treatment), Time = df2 |> 
      dplyr::filter(Time == "3h") |> 
      pull(Time),
    col = ann_colors,
    show_annotation_name = F,
    show_legend = c("Treatment" = F, "Time" = F))
  
  ht2<-ComplexHeatmap::Heatmap(
    scaled_mat[,colnames(scaled_mat)[grepl("t2", colnames(scaled_mat))]],
    col = col_fun_brgn,
    column_title = "            Samples",
    column_title_side = "bottom",
    row_title = "Genes",
    show_row_names = F,
    show_column_names = F,
    column_names_gp = gpar(fontsize = 8),
    #column_labels = newLabs$NewName,
    row_names_gp = gpar(fontsize = 4),
    show_row_dend = F,
    cluster_columns = F,
    cluster_rows = T,
    heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")),
    top_annotation = ha2, 
    row_km_repeats = 100,
    column_km_repeats = 100,
    show_heatmap_legend = F
    #column_split = 2
  )
  
  draw(ht2)
  
  ha3<-HeatmapAnnotation(
    Treatment = df2 |> 
      dplyr::filter(Time == "8h") |> 
      pull(Treatment), Time = df2 |> 
      dplyr::filter(Time == "8h") |> 
      pull(Time),
    col = ann_colors,
    show_annotation_name = F,
    show_legend = c("Treatment" = F, "Time" = F))
  
  ht3<-ComplexHeatmap::Heatmap(
    scaled_mat[,colnames(scaled_mat)[grepl("t3", colnames(scaled_mat))]],
    col = col_fun_brgn,
    column_title = "",
    column_title_side = "bottom",
    row_title = "",
    show_row_names = F,
    show_column_names = F,
    column_names_gp = gpar(fontsize = 8),
    #column_labels = newLabs$NewName,
    row_names_gp = gpar(fontsize = 4),
    show_row_dend = F,
    cluster_columns = F,
    cluster_rows = T,
    heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")),
    top_annotation = ha3, 
    row_km_repeats = 100,
    column_km_repeats = 100,
    show_heatmap_legend = F
    #column_split = 2
  )
  
  draw(ht3)
  #ht4
  ha4<-HeatmapAnnotation(
    Treatment = df2 |> 
      dplyr::filter(Time == "24h") |> 
      pull(Treatment), Time = df2 |> 
      dplyr::filter(Time == "24h") |> 
      pull(Time),
    col = ann_colors,
    annotation_legend_param = list(Treatment = list(title = "Treatment"), 
                                   Time = list(title = "Time", at = times)),
    show_annotation_name = F)
  
  ht4<-ComplexHeatmap::Heatmap(
    scaled_mat[,colnames(scaled_mat)[grepl("t4", colnames(scaled_mat))]],
    col = col_fun_brgn,
    column_title = "",
    column_title_side = "bottom",
    row_title = "Genes",
    show_row_names = F,
    show_column_names = F,
    column_names_gp = gpar(fontsize = 8),
    #column_labels = newLabs$NewName,
    row_names_gp = gpar(fontsize = 4),
    show_row_dend = F,
    cluster_columns = F,
    cluster_rows = T,
    heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")),
    top_annotation = ha4, 
    row_km_repeats = 100,
    column_km_repeats = 100,
    #column_split = 2
  )
  
  draw(ht4)
  ht_list = ht1 + ht2 + ht3 + ht4
  draw(ht_list)
  
# save
  filepath<-paste0("publication/SOM-Fig4-", sampSet)
  svglite::svglite(paste0(filepath, ".svg"),width = 9, height = 7)
  draw(ht_list)
  
dev.off()
}else if(sampSet=="29samp"){
  # Order time points correctly for plots
  lfc<-2
  padj<-0.05
  times <- c("3h", "8h", "24h")
  
  ann_colors = list(
    Treatment = c(DMSO = trt_vir_cols[1], HHQ = trt_vir_cols[2]),
    Time = c(
      `3h` = grays[3],
      `8h` = grays[4],
      `24h` = grays[5]
    )
  )
  
  vir<-readRDS(paste0(res_path, "deseq_glm_CountsContsOfInterest_", sampSet, "_virus.rds"))
  virDegs1<-lapply(vir, function(x) x |> dplyr::filter(padjIHW<=padj))
  unlist(lapply(virDegs1, function(x) nrow(x)))
  virDegs2<-lapply(virDegs1, function(x) x |> dplyr::filter(abs(log2FCshrink_ashr)>=lfc))
  unlist(lapply(virDegs2, function(x) nrow(x)))
  candidate_genes<-unique(unname(unlist(lapply(virDegs, function(x) x$locus_tag))))
  sample_info<-read_tsv(paste0(base_path, "data/coldata_", sampSet, ".txt"))
  virAnnot<-readxl::read_xlsx(paste0(base_path, "data/Annotations-virus/annot-virus-ext.xlsx"), sheet = "All")
  candidate_genes2<-virAnnot |> dplyr::filter(locus_tag %in% candidate_genes) |> pull(accession)
  
# MyPheatmapClust(
#   mat,
#   main = "Test",
#   scale = "row",
#   cluster_rows = T,
#   show_rownames = F,
#   cluster_cols = F,
#   color = htMapCols,
#   annotation_col = df2[,-3],
#   gaps_col = c(9,19),
#   annotation_colors = ann_colors,
#   cutree_rows = 5
# )
pheatmap::pheatmap(scaled_mat,
                    main = "Test",
                    cluster_rows = T,
                    show_rownames = F,
                    cluster_cols = F,
                    color = br_gn,
                    annotation_col = df2[,-3],
                    gaps_col = c(9,19), annotation_colors = ann_colors)


  
# plot complex heatmaps
#ha<-HeatmapAnnotation(
#  Treatment = df2$Treatment, Time = df2$Time,
#  col = ann_colors,
#  annotation_legend_param = list(Treatment = list(title = "Treatment"), 
#                                 Time = list(title = "Time", at = times)),
#  annotation_name_side = "left",which = "column")

ht<-ComplexHeatmap::Heatmap(
  scaled_mat,
  col = col_fun_brgn,
  column_title = "Samples",
  column_title_side = "bottom",
  row_title = "Genes",
  show_row_names = F,
  show_column_names = F,
  column_names_gp = gpar(fontsize = 8),
  #column_labels = newLabs$NewName,
  row_names_gp = gpar(fontsize = 4),
  show_row_dend = F,
  cluster_columns = F,
  cluster_rows = T,
  heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")),
  top_annotation = ha, 
  row_km_repeats = 100,
  column_km_repeats = 100,
  split = 3,gap = unit(3, "mm"),
  show_heatmap_legend = T
  #column_split = 2
)
draw(ht)
htb<-ComplexHeatmap::Heatmap(
  scaled_mat,
  col = col_fun_brgn,
  column_title = "Samples",
  column_title_side = "bottom",
  row_title = "Genes",
  show_row_names = F,
  show_column_names = F,
  column_names_gp = gpar(fontsize = 8),
  #column_labels = newLabs$NewName,
  row_names_gp = gpar(fontsize = 4),
  show_row_dend = F,
  cluster_columns = T,
  cluster_rows = T,
  heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")),
  top_annotation = ha, 
  row_km_repeats = 100,
  column_km_repeats = 100,
  split = 3,gap = unit(3, "mm"),
  column_split = 3,
  show_heatmap_legend = T
  #column_split = 2
)

ha1<-HeatmapAnnotation(
  Treatment =df2 |> dplyr::filter(Time == "3h") |> pull(Treatment), Time = df2 |> dplyr::filter(Time == "3h") |> pull(Time),
  col = ann_colors,
 # annotation_legend_param = 
#    list(Treatment = list(title = "Treatment"), 
#    Time = list(title = "Time", at = times)),
  annotation_name_side = "left",which = "column",
show_legend = c("Treatment" = F, "Time" = F))

ht1<-ComplexHeatmap::Heatmap(
  scaled_mat[,colnames(scaled_mat)[grepl("t2", colnames(scaled_mat))]],
  col = col_fun_brgn,
  column_title = "",
  column_title_side = "bottom",
  row_title = "Genes",
  show_row_names = F,
  show_column_names = F,
  column_names_gp = gpar(fontsize = 8),
  #column_labels = newLabs$NewName,
  row_names_gp = gpar(fontsize = 4),
  show_row_dend = T,
  cluster_columns = F,
  cluster_rows = T,
  heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")),
  top_annotation = ha1, 
  row_km_repeats = 100,
  column_km_repeats = 100,
  split = 3,
  gap = unit(3, "mm"),
  show_heatmap_legend = F
  #column_split = 2
  )

draw(ht1)

ha2<-HeatmapAnnotation(
  Treatment =df2 |> 
    dplyr::filter(Time == "8h") |> 
    pull(Treatment), Time = df2 |> 
    dplyr::filter(Time == "8h") |> 
    pull(Time),
  col = ann_colors,
 #annotation_legend_param = list(Treatment = list(title = "Treatment"), 
#                                 Time = list(title = "Time", at = times)),
  #annotation_name_side = "left",which = "column",
show_annotation_name = F,
show_legend = c("Treatment" = F, "Time" = F))

ht2<-ComplexHeatmap::Heatmap(
  scaled_mat[,colnames(scaled_mat)[grepl("t3", colnames(scaled_mat))]],
  col = col_fun_brgn,
  column_title = "Samples",
  column_title_side = "bottom",
  row_title = "Genes",
  show_row_names = F,
  show_column_names = F,
  column_names_gp = gpar(fontsize = 8),
  #column_labels = newLabs$NewName,
  row_names_gp = gpar(fontsize = 4),
  show_row_dend = F,
  cluster_columns = F,
  cluster_rows = T,
  heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")),
  top_annotation = ha2, 
  row_km_repeats = 100,
  column_km_repeats = 100,
  show_heatmap_legend = F
  #column_split = 2
)

draw(ht2)
ha3<-HeatmapAnnotation(
  Treatment = df2 |> 
    dplyr::filter(Time == "24h") |> 
    pull(Treatment), Time = df2 |> 
    dplyr::filter(Time == "24h") |> 
    pull(Time),
  col = ann_colors,
  annotation_legend_param = list(Treatment = list(title = "Treatment"), 
                                 Time = list(title = "Time", at = times)),
  show_annotation_name = F)

ht3<-ComplexHeatmap::Heatmap(
  scaled_mat[,colnames(scaled_mat)[grepl("t4", colnames(scaled_mat))]],
  col = col_fun_brgn,
  column_title = "",
  column_title_side = "bottom",
  row_title = "Genes",
  show_row_names = F,
  show_column_names = F,
  column_names_gp = gpar(fontsize = 8),
  #column_labels = newLabs$NewName,
  row_names_gp = gpar(fontsize = 4),
  show_row_dend = F,
  cluster_columns = F,
  cluster_rows = T,
  heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")),
  top_annotation = ha3, 
  row_km_repeats = 100,
  column_km_repeats = 100,
  #column_split = 2
)

draw(ht3)
ht_list = ht1+ ht2+ht3
draw(ht_list)
#res = 150L

filepath<-paste0(res_path, "figs/VirZscoreHeatMap_", sampSet)
svglite::svglite(paste0(filepath, ".svg"),width = 9, height = 7)
draw(ht_list)
dev.off()
}

# Render the svg into a png image with rsvg via magick
res=150
img <- magick::image_read_svg(paste0(filepath, ".svg"), width = 1080)
png(paste0(filepath, ".png"),
    width = 1200,
    height = 1000,
    res = res)
grid::grid.newpage()
draw(ht_list)
dev.off()

# SOM Fig S5 ---------
mat <- assay(rld_blind)
metadata <- colData(dds)[,-1]
library(PCAtools)
nprop <-0.95
ntop <- round(nrow(dds)*nprop)

annot<-annot<-readRDS("data/Annotations-virus/annot-virus-edgeR.rds") |>
  dplyr::select(tx_id, locus_tag, rowLabel, everything())

p <- PCAtools::pca(
  mat = mat,
  metadata = metadata,
  removeVar = 1-nprop) # fraction of low variance data to exclude
df2 <- as.data.frame(colData(dds)[,c("treatment","timePt")])


ld <- as_tibble(getLoadings(p), rownames = "#locus_tag")



rownames(p$loadings) <- as_tibble(rownames(p$loadings)) |>
  setNames("rowLabel") |>
  left_join(annot) |>
  dplyr::select(c("rowLabel")) |>
  deframe() 

times <- c("45min", "3h", "8h", "24h")

pca_data12<-DESeq2::plotPCA(
  rld_blind, ntop =  ntop, 
  #vsd_blind, ntop =  ntop, 
  returnData = T, intgroup = c("trt_vir", "timePt"), pcsToUse=1:2)

pca_data23<-DESeq2::plotPCA(
  rld_blind, ntop =  ntop,
  #vsd_blind, ntop =  ntop, 
  returnData = T, intgroup = c("trt_vir", "timePt"), pcsToUse=c(2,3))


#fix names of treatments and time points for plots
percentVar12 <- round(100 * attr(pca_data12, "percentVar"))
percentVar23 <- round(100 * attr(pca_data23, "percentVar"))

(pca_data12$trt_vir<-PrepNamesForPlots(trt_vir = pca_data12$trt_vir)$trt_vir)
(pca_data23$trt_vir<-PrepNamesForPlots(trt_vir = pca_data23$trt_vir)$trt_vir)
#fix time points
(pca_data12$timePt2<-PrepNamesForPlots(timePt = pca_data12$timePt)$timePt)
(pca_data23$timePt2<-PrepNamesForPlots(timePt = pca_data23$timePt)$timePt)


pca_data12$timePt2<-factor(pca_data12$timePt2, levels = times, ordered = T)
pca_data23$timePt2<-factor(pca_data23$timePt2, levels = times, ordered = T)


trt_vir_cols<-c(paired_cols[1:2], paired_cols[7:8])
trt_vir_cols<- trt_vir_cols[c(2,4)]

timept_shapes <- c(19, 15, 18, 17)
size_vals <- c(3, 3, 4, 3)
size_vals2<-c(3, 3, 4, 3)

pplot1a <- pca_data12 |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(
    shape = timePt2,
    colour = trt_vir,
    size = timePt2),
    alpha = 0.8
  ) +
  #make it pretty
  scale_color_manual(values = trt_vir_cols) +
  scale_shape_manual(values = timept_shapes) +
  coord_fixed() +
  #make all points same size
  scale_size_manual(values = size_vals) +
  #ggtitle("PCA of Variance Stabilized Counts" )+
  xlab(paste0("PC1 (", percentVar12[1], "% variance)")) +
  ylab(paste0("PC2 (", percentVar12[2], "% variance)")) +
  bb_theme() +
  guides(shape = guide_legend(override.aes = list(size = size_vals2))) +
  guides(colour = guide_legend(override.aes = list(size = length(size_vals))))
#pplot1a
pplot1a_nolg <- pplot1a + theme(legend.position = "none")

#pplot1a

pplot1b <- pca_data23 |>
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(
    shape = timePt2,
    colour = trt_vir,
    size = timePt2),
    alpha=0.8
  ) +
  #make it pretty
  scale_color_manual(values = trt_vir_cols) +
  scale_shape_manual(values = timept_shapes) +
  coord_fixed() +
  #make all points same size
  scale_size_manual(values = size_vals) +
  #ggtitle("PCA of Variance Stabilized Counts" )+
  xlab(paste0("PC2 (", percentVar23[1], "% variance)")) +
  ylab(paste0("PC3 (", percentVar23[2], "% variance)")) +
  bb_theme() +
  guides(shape = guide_legend(override.aes = list(size = size_vals2))) +
  guides(colour = guide_legend(override.aes = list(size = length(size_vals))))

#prettier PCAs
p1 <- pplot1a_nolg / pplot1b + plot_annotation(title = "Pretty PCA, no ellipses", theme =
                                                 bb_theme())
p1
#part C
salmon_gene_quants<-read_tsv(paste0(res_path, "tables/salmon-output-processed-filt", ext,".tsv.gz")) 
tib <- tibble(sample_name = unique(salmon_gene_quants$sample_name),SplitSampleName(unique(salmon_gene_quants$sample_name)))
genes_keep_final<-readRDS(paste0(res_path, "genes_keep_final2", ext, ".rds"))
salmon_gene_quants_filt <- 
  salmon_gene_quants |>
  dplyr::filter(tx_id %in% (
    annot |>
      dplyr::filter(locus_tag %in% genes_keep_final) |>
      pull(tx_id)
  ))
salmon_gene_quants_filt <- left_join(salmon_gene_quants_filt, tib)

partC<-salmon_gene_quants_filt |>
  group_by(timePt) |> 
  summarise(`Reads Mapped (Millions)` = round(sum(n_reads),2)) |> 
  ungroup() |> 
  dplyr::rename(`Time Point` =  timePt) |>
  mutate(`Time Point` = c("45 min", "3 hr", "8 hr", "24 hr")) |>
  mutate(Prop = round(`Reads Mapped (Millions)`/sum(`Reads Mapped (Millions)`),7)) |>
  mutate(`% of total` = round(Prop*100, 2)) |>
  select(-Prop) |>
  gt() |>  opt_stylize(add_row_striping = TRUE, style = 1) |>
  opt_table_font(
    font = list(
      google_font(name = "Roboto"),
      "Cochin", "serif"
    )) |>
  fmt_scientific(columns = `Reads Mapped (Millions)`)
library(patchwork)

filepath<-paste0("publication/SOM-Fig5-PartC-", sampSet)
svglite::svglite(paste0(filepath, ".svg"),width = 9, height = 7)
salmon_gene_quants_filt |>
  group_by(timePt) |> 
  summarise(`Reads Mapped (Millions)` = round(sum(n_reads),2)) |> 
  ungroup() |> 
  dplyr::rename(`Time Point` =  timePt) |>
  mutate(`Time Point` = c("45 min", "3 hr", "8 hr", "24 hr")) |>
  mutate(Prop = round(`Reads Mapped (Millions)`/sum(`Reads Mapped (Millions)`),7)) |>
  mutate(`% of total` = round(Prop*100, 2)) |>
  select(-Prop) |>
  gt() |>  opt_stylize(add_row_striping = TRUE, style = 1) |>
  opt_table_font(
    font = list(
      google_font(name = "Roboto"),
      "Cochin", "serif"
    )) |>
  fmt_scientific(columns = `Reads Mapped (Millions)`)
dev.off()
#C
partC |>
  gtsave(paste0("publication/SOM-Fig5-C-", sampSet, ".png"), expand = 10)
svglite::svglite(paste0(filepath, ".svg"),width = 9, height = 7)
partC
dev.off()
#A-B
filepath<-paste0("publication/SOM-Fig5-AB-", sampSet)
svglite::svglite(paste0(filepath, ".svg"),width = 9, height = 7)
p1
dev.off()

res=150
img <- magick::image_read_svg(paste0(filepath, ".svg"), width = 1080)
png(paste0(filepath, ".png"),
    width = 1200,
    height = 1000,
    res = res)
grid::grid.newpage()
p1
dev.off()
# save

# SOM Data File 5 ---------


annot<-readRDS("data/Annotations-host/annot-host-ext2-2026-06-28.rds")
annot<-annot |> distinct()
annot<-annot |> select(-c(EnsemblProtists, origin))
annot<-annot |> distinct()

legend<-tibble(`Column Name` = colnames(host_annot), 
               `Description` = c("Genkank transcript ID","locus tag (ncbi). Note, if more than one transcript exists for a given locus tag, a '-' was added at the end followed by a number to differentiate them"))

# Playground ----------
MyPheatmapClust(
  scaled_mat,
  main = "Test",
  cluster_rows = T,
  show_rownames = F,
  cluster_cols = F,
  color = br_gn,
  annotation_col = df2[,-3],
  cutree_rows = 5
)


#https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html

# Som Fig S4
# Summarise counts 
trans_cts_mean <- trans_cts %>% 
  # convert to long format
  pivot_longer(cols = DMSO_inf_t2:HHQ_inf_t4_5, names_to = "names", values_to = "cts")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("names")) %>% 
  # filter to retain only genes of interest
  filter(gene %in% candidate_genes2) %>% 
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(gene,treatment, timePt) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  na.omit() |>
  ungroup()

trans_cts_mean <- trans_cts_mean |>
  mutate(timePt = PrepNamesForPlots(timePt=timePt)$timePt)

hclust_matrix <- trans_cts[,-c(1:17)] %>% 
  dplyr::select(-"gene") |>
  as.matrix()

# assign rownames
rownames(hclust_matrix) <- trans_cts$gene
hclust_matrix <- hclust_matrix[candidate_genes2, ]

hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")
plot(gene_hclust)
abline(h = 3, col = "brown", lwd = 2) #
abline(h = 5, col = "brown", lwd = 2) #

gene_cluster <- cutree(gene_hclust, k = 3) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(gene = name, cluster = value)

head(gene_cluster)

trans_cts_cluster <- trans_cts_mean %>% 
  inner_join(gene_cluster, by = "gene")

head(trans_cts_cluster)

library(viridis)
library(colorspace)
library(cols4all)
virid<-c4a(palette="viridis")
ziss<c4a(palette="Zissou 1")
names(ziss)<-c(1:length(ziss))
trans_cts_cluster<-trans_cts_cluster %>% 
  mutate(time = factor(parse_number(timePt), ordered = T))|>
  mutate(cluster = factor(cluster, ordered = T))

p<-trans_cts_cluster |>   ggplot(aes(time, mean_cts_scaled, label = gene, group=cluster)) +
  geom_line(aes(group = gene, color = cluster, alpha = 0.7), linewidth = 1, show.legend = F) +
  facet_grid(rows = vars(cluster), cols = vars(treatment)) +
  xlab("Hours") + 
  ylab("Scaled Mean Counts") +
  #ggrepel::geom_text_repel(aes())+
  scale_colour_manual(values = virid) +
  #viridis::scale_color_viridis(discrete = TRUE)+
  bb_theme() 

degsHT<-Heatmap(hclust_matrix, show_row_names = FALSE, col = col_fun_brgn, cluster_columns = F, split = 3, heatmap_legend_param = list(title = "", at = c(-3,-2, -1, 0, 1,2,3),legend_height = unit(6, "cm")))
draw(degsHT)


# wgnca
#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#44_Perform_DESeq2_normalization_and_transformation




# Figure 3 ------------
subs <- "host"
sampSet <- "76samp" #
res_path <- paste0(base_path, "host/76samp/")
(ext <- paste0("_", sampSet, "_", subs))
#timeP4upsetDeseqINTERSECT.png
#timeP3upsetDeseqINTERSECT.png
#timeP2upsetDeseqINTERSECT.png
#timePt1upsetDeseqINTERSECT.png
res<-readRDS(paste0(res_path, "deseq_glm_CountsContsOfInterest", ext, ".rds"))
l.names <- names(res)
gc()
l.t1names <- l.names[grepl("t1", l.names)]
l.t2names <- l.names[grepl("t2", l.names)]
l.t3names <- l.names[grepl("t3", l.names)]
l.t4names <- l.names[grepl("t4", l.names)]

res2 <- vector('list', length(res))
names(res2) <- names(res)
res2 <- mclapply2(res, function(x) {
  x |>
    mutate(
      Ranking = ifelse(log2FoldChange < 0, -1, ifelse(log2FoldChange > 0, 1, 0)) *
        -log10(pvalueRaw),
      .after = locus_tag
    ) |>
    mutate(DE = padjIHW <= padj&
             abs(log2FCshrink_ashr) >= lfc,
           .after = Ranking) |>
    group_by(locus_tag) |>
    dplyr::slice(1) |>
    arrange(desc(Ranking))
})

res2t1 <- res2[l.t1names]
res2t2 <- res2[l.t2names]
res2t3 <- res2[l.t3names]
res2t4 <- res2[l.t4names]
temp <- as_tibble(do.call(rbind, lapply(names(res2t4), function(x)
  str_split_1(x, "v")))) |> mutate(V1 = PrettyTrtVir(V1), V2 = PrettyTrtVir(V2)) |>
  mutate(contrast = paste0(V1, "\nx\n", V2))
names(res2t4) <- temp$contrast
names(res2t3) <- temp$contrast
names(res2t2) <- temp$contrast
names(res2t1) <- temp$contrast
res2t4 <- lapply(res2t4, function(x)
  x |> dplyr::filter(DE == T) |> pull(locus_tag))
res2t3 <- lapply(res2t3, function(x)
  x |> dplyr::filter(DE == T) |> pull(locus_tag))
res2t2 <- lapply(res2t2, function(x)
  x |> dplyr::filter(DE == T) |> pull(locus_tag))
res2t1 <- lapply(res2t1, function(x)
  x |> dplyr::filter(DE == T) |> pull(locus_tag))


filepath <- paste0(res_path, "figs/timeP", 1:4, "upsetDeseqINTERSECT.svg")
titles <- c("45min", "3h", "8h", "24h")
full_list <- list(t1 = res2t1,
                  t2 = res2t2,
                  t3 = res2t3,
                  t4 = res2t4)
#plot
res = 300L
for (t in 1:4) {
  svglite::svglite(filepath[t], width = 9, height = 7)
  MyUpsetPlot(
    x = full_list[[t]],
    mode = "intersect",
    main = titles[t],
    file = filepath[t] ,
    png = F
  )
  dev.off()
}
# Render the svg into a png image with rsvg via magick
img <- magick::image_read_svg(paste0(filepath, ".svg"), width = 1080)
magick::image_write(img, paste0(filepath, ".png"))

