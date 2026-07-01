# Plots for Publication
source("dge/scripts/Functions.R")
source("dge/scripts/plot_funcs.R")
library(conflicted)
conflicts_prefer(dplyr::filter)
## this font is finicky
#systemfonts::fonts_as_import(family = "Roboto Condensed")
#systemfonts::match_font("Roboto Condensed")
#base_path <- "~/Library/CloudStorage/GoogleDrive-barbarabitarello@gmail.com/My Drive/BitaLab_not_shared/Research/rna-seq-host-virus/data_and_res_gh_repo/"
base_path<-"~/Documents/GitHub/Ehux_Ehv201/scratch/"


# Template ----
if(FALSE){ #not run
plot_path = getwd()
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


# Figure 3 ------------
subs <- "host"
sampSet <- "76samp" #
res_path <- paste0(base_path, "host/76samp/")
(ext <- paste0("_", sampSet, "_", subs))
#timeP4upsetDeseqINTERSECT.png
#timeP3upsetDeseqINTERSECT.png
#timeP2upsetDeseqINTERSECT.png
#timePt1upsetDeseqINTERSECT.png

res <- readRDS(paste0(res_path, "deseq.tb.list.contsOfInterest.rds"))
res <- res[-c(1:4)]
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
    mutate(DE = padjIHW <= 0.05 &
             abs(log2FCshrink_ashr) >= 1,
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


# SOM S2: Heatplots for virus ------

# 12/26: replotting with all time points per KW's request by email.

subs <- "virus"
sampSet <- "39samp"
#sampSet <- "29samp"

# outfile path
res_path <- paste0(base_path, subs, "/", sampSet, "/")
(ext <- paste0("_", sampSet, "_", subs))

# read in data
dds <- readRDS(paste0(res_path, "deseq_dds", ext, ".rds"))
rld_blind <- readRDS(paste0(res_path, "rld_blind", ext, ".rds"))
vsd_blind <- readRDS(paste0(res_path, "vsd_blind", ext, ".rds"))
rld <- readRDS(paste0(res_path, "rld", ext, ".rds"))
vsd <- readRDS(paste0(res_path, "vsd", ext, ".rds"))

#colors
trt_vir_cols <- c(paired_cols[2], paired_cols[8])

#if only t2-t4:
if(sampSet=="29samp"){
# Order time points correctly for plots
times <- c("3h", "8h", "24h")

ann_colors = list(
  Treatment = c(DMSO = trt_vir_cols[1], HHQ = trt_vir_cols[2]),
  Time = c(
    `3h` = grays[3],
    `8h` = grays[4],
    `24h` = grays[5]
  )
)
}else if(sampSet=="39samp"){
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
}

# parse data
test_result<-readRDS(paste0(res_path, "deseq_glm_CountsContsOfInterest_", sampSet, "_virus.rds"))
test_result<-lapply(test_result, function(x) x |> dplyr::filter(padjBH<=0.05))
candidate_genes<-unique(unname(unlist(lapply(test_result, function(x) x$locus_tag))))
sample_info<-read_tsv(paste0(base_path, "data/coldata_", sampSet, ".txt"))
virAnnot<-readxl::read_xlsx(paste0(base_path, "data/Annotations-virus/annot-virus-ext.xlsx"), sheet = "All")
candidate_genes2<-virAnnot |> dplyr::filter(locus_tag %in% candidate_genes) |> pull(accession)
trans_cts<-readxl::read_xlsx(paste0(res_path, "tables/deseq_counts_contrastsOfInterest_", sampSet, "_virus.xlsx"), sheet = "variance_stabilized_counts")
trans_cts<-trans_cts |> mutate(gene = accession)
#trans_cts |> filter(locus_tag %in% keep_genes_final)
library(DESeq2)
#mat<- assay(rld_blind)
mat <- assay(rld)
#mat <- assay(vsd_blind)
rownames(mat)<-trans_cts[match(rownames(mat),trans_cts$locus_tag),]$gene
scaled_mat = t(scale(t(mat)))

coldata <- read_tsv(paste0(base_path, "data/coldata_", sampSet, ".txt"),
                    show_col_types = FALSE)
df2 <- as.data.frame(colData(dds)[, c("treatment", "timePt")])
colnames(df2)[1] <- 'trt'
df2$timePt <- PrepNamesForPlots(timePt = df2$timePt)$timePt
# arrange by time
df2$timePt2 <- parse_number(gsub("h", "", df2$timePt))
df2<-df2 |> mutate(timePt2=ifelse(timePt2==45, 0.75, timePt2))
df2 <- df2 |> arrange(timePt2) |> dplyr::select(-timePt2)
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

if(sampSet=="29samp"){
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
}else{

  pheatmap::pheatmap(scaled_mat,
                     main = "Test",
                     cluster_rows = T,
                     show_rownames = F,
                     cluster_cols = F,
                     color = br_gn,
                     annotation_col = df2[,-3],
                     gaps_col = c(10,19,29), annotation_colors = ann_colors)
}
#check
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

# plot complex heatmaps
#ha<-HeatmapAnnotation(
#  Treatment = df2$Treatment, Time = df2$Time,
#  col = ann_colors,
#  annotation_legend_param = list(Treatment = list(title = "Treatment"), 
#                                 Time = list(title = "Time", at = times)),
#  annotation_name_side = "left",which = "column")
if(sampSet=="29samp"){
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
res = 150L

# save
filepath<-paste0(res_path, "figs/VirZscoreHeatMap")
svglite::svglite(paste0(filepath, ".svg"),width = 9, height = 7)
draw(ht_list)
dev.off()
}else{
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
res = 150L

# save
filepath<-paste0(res_path, "figs/VirZscoreHeatMap_", sampSet)
svglite::svglite(paste0(filepath, ".svg"),width = 9, height = 7)
draw(ht_list)
dev.off()
}
# Render the svg into a png image with rsvg via magick
img <- magick::image_read_svg(paste0(filepath, ".svg"), width = 1080)
png(paste0(filepath, ".png"),
    width = 1200,
    height = 1000,
    res = res)
grid::grid.newpage()
draw(ht_list)
dev.off()

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
