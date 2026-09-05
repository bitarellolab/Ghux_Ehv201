#Functions

##remove NA columns

NAColOmit<-function(x = NULL){
  #x is a tibble/data.frame/data.table
  require(tidyverse)
  x<- as_tibble(x)
  n1 <- nrow(x)
  j1 <- ncol(x)
  rem_these <- names(which(apply(x, 2, function(x) sum(is.na(x)) == n1)))
  paste0("removed cols: ", paste0(rem_these, collapse = ", "))
  
  x |> dplyr::select(-all_of(rem_these))
}

NAPerRow<-function(x = NULL){
  #x is a tibble/data.frame/data.table
  require(tidyverse)
  x<- as_tibble(x)
  n1 <- nrow(x)
  j1 <- ncol(x)
  #rem_these <- names(which(apply(x, 2, function(x) sum(is.na(x)) == n1)))
  totalNAs<-apply(x, 1, function(x) sum(is.na(x)))
  x |> dplyr::mutate(TotalNAs = totalNAs)
}
# for dplyr
allNAs<-function(x){sum(is.na(x))}

## meanSDPlot (modified):
# modified from this version https://gist.github.com/mikelove/0a3c19a8512fb36453078df48f20d283
MyMeanSdPlot <- function(x,
                         ranks = TRUE,
                         xlab = ifelse(ranks, "rank(mean)", "mean"),
                         ylab = "sd",
                         pch,
                         plot = TRUE,
                         bins = 50,
                         main = "GIVE ME A TITLE",
                         ...) {
  stopifnot(is.logical(ranks), length(ranks) == 1, !is.na(ranks))
  
  n = nrow(x)
  if (n == 0L) {
    warning("In 'meanSdPlot': input matrix 'x' has 0 rows. There is nothing to be done.")
    return()
  }
  if (!missing(pch)) {
    warning("In 'meanSdPlot': 'pch' is ignored.")
  }
  
  px   = rowMeans(x, na.rm = TRUE)
  ## here I (Mike Love) edit rowV to rowVars
  py   = sqrt(matrixStats::rowVars(x, mean = px, na.rm = TRUE))
  rpx  = rank(px, na.last = FALSE, ties.method = "random")
  
  ## run median with centers at dm, 2*dm, 3*dm,... and width 2*dm
  dm        = 0.025
  midpoints = seq(dm, 1 - dm, by = dm)
  within    = function(x, x1, x2) {
    (x >= x1) & (x <= x2)
  }
  mediwind  = function(mp)
    median(py[within(rpx / n, mp - 2 * dm, mp + 2 * dm)], na.rm = TRUE)
  rq.sds    = sapply(midpoints, mediwind)
  
  res = if (ranks) {
    list(
      rank = midpoints * n,
      sd = rq.sds,
      px = rpx,
      py = py
    )
  } else {
    list(
      quantile = quantile(px, probs = midpoints, na.rm = TRUE),
      sd = rq.sds,
      px = px,
      py = py
    )
  }
  
  fmt = function()
    function(x)
      format(round(x, 0), nsmall = 0L, scientific = FALSE)
  
  res$gg = ggplot(data.frame(px = res$px, py = res$py),
                  aes_string(x = "px", y = "py")) + xlab(xlab) + ylab(ylab) +
    geom_hex(bins = bins, ...) +
    scale_fill_gradient(name = "count",
                        trans = "log",
                        labels = fmt()) +
    geom_line(
      aes_string(x = "x", y = "y"),
      data = data.frame(x = res[[1]], y = res$sd),
      color = "red"
    ) +
    ggtitle(label = main) + theme_minimal()
  
  if (plot)
    print(res$gg)
  
  return(invisible(res))
}

#extract variables from sample name ----

SplitSampleName <- function(sample = "DMSO_cntl_t3_2") {
  time <- str_extract(sample, "[\\d]")
  trt <- str_extract(sample, "^[^_]+")
  virus <- str_extract(sample, "(cntl|inf)")
  rep <- str_extract(sample, "[\\d]$")
  res <- tibble(
    timePt = time,
    trt = trt,
    virus = virus,
    rep = rep
  )
  return(res)
}

#extract info from contrast name -----
 SplitContrastName <-
   function(contrast = "HHQ_inf_t1vDMSO_cntl_t3", sep = "v") {
     l<-length(contrast)
     conts <- data.table::tstrsplit(contrast, sep, fixed = T)
     if(l == 1){
     names(conts) <- c("expL", "cntlL")
     conts2 <- lapply(conts, function(x)
       SplitSampleName(sample = x))
     names(conts2) <- c("expL", "cntlL")     
     res<-do.call(rbind, lapply(1:2, function(x) conts2[[x]] |> mutate(contrastL = names(conts)[x], contrast = contrast)))
     res <- do.call(rbind, lapply(1:2, function(x)
       conts2[[x]])) |> mutate(contrast = contrast)
     }else{
     for(i in 1:l){
       names(conts[[i]])<-c("expL", "cntlL")
     }
     lapply(conts, function(x) lapply(x, function(y) separate()))
     }
     return(res)
   }
#


# make trt_vir and timePt names ready for plots ----
PrepNamesForPlots <- function(trt_vir = NULL, timePt = NULL) {
  time_1 <- gsub("t4", "24h", gsub("t3", "8h", gsub("t2", "3h", gsub(
    "t1", "45min", timePt
  ))))
  trt_vir_1 <- paste0(gsub("_", " (", trt_vir), ")")
  return(list(timePt = time_1, trt_vir = trt_vir_1))
}



# remove reps form sample name ----
PrepNamesForPlots2 <- function(names = "DMSO_inf_t1_1") {
  name_1 <- gsub("_(1|2|3|4|5)", "", names)
  return(name_1)
}

#remove reps and inf/cntlfrom sample name (for virus)----
PrepNamesForPlots3 <- function(names = "DMSO_inf_t1_1") {
  name_1 <- gsub("_(cntl|inf)", "", gsub("_(1|2|3|4|5)", "", names))
  return(name_1)
}


#make better labels

PrettyTrtVir <- function(x) {
  tibble(orig = x)  |> separate(orig, c("trt", "virus"), "_", remove = F) |> mutate(new = paste0(trt, " (", virus, ")")) |> dplyr::select(new) |> deframe()
  
}

# custom version of parallel::mclapply that shoes progress bar. ----
mclapply2 <- function(X,
                      FUN,
                      ...,
                      mc.preschedule = FALSE,
                      mc.set.seed = TRUE,
                      mc.silent = FALSE,
                      mc.cores = getOption("mc.cores", 2L),
                      mc.cleanup = TRUE,
                      mc.allow.recursive = TRUE,
                      mc.progress = TRUE,
                      mc.style = 3)
{
  require(parallel)
  if (!is.vector(X) || is.object(X))
    X <- as.list(X)
  
  if (mc.progress) {
    f <- fifo(tempfile(), open = "w+b", blocking = T)
    p <- parallel:::mcfork()
    pb <- txtProgressBar(0, length(X), style = mc.style)
    setTxtProgressBar(pb, 0)
    progress <- 0
    if (inherits(p, "masterProcess")) {
      while (progress < length(X)) {
        readBin(f, "double")
        progress <- progress + 1
        setTxtProgressBar(pb, progress)
      }
      cat("\n")
      parallel:::mcexit()
    }
  }
  tryCatch({
    result <- mclapply(X, ..., function(...) {
      res <- FUN(...)
      if (mc.progress)
        writeBin(1, f)
      res
    }, mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.silent = mc.silent, mc.cores = mc.cores, mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive)
    
  }, finally = {
    if (mc.progress)
      close(f)
  })
  result
}


## check if any GO annot

isGOannot <- function(x = x,
                      cols = c("GO_BP", "GO_CC", "GO_MF", "GO_GO", "GO_IDs")) {
  x |> mutate(GOAnnotAvail = rowSums(is.na(x |> dplyr::select({
    cols
  })) == F, na.rm = T) != 0)
}

#isGOannot <- function(x = x,
#                      cols = c("GO_BP", "GO_CC", "GO_MF", "GO_GO", "GO_IDs", "GOs")) {
#  x |> mutate(GOAnnotAvail = rowSums(is.na(x |> dplyr::select({
#    cols
#  })) == F, na.rm = T) != 0)
#}


# Find Top Degs and plot Volcano ----
VolcanoFunc <- function(cont = NULL,
                        padj = 0.05,
                        lfc = c(2, 1.5), subs = "virus") {
  message("contrast: ", cont)
 # message("lfc: ", lfc)
  df <- glm_subset[[cont]] |>
    filter(contrast == cont) |>
    arrange(desc(abs(log2FCshrink_ashr))) |>
    mutate(thresh_padjlfc = padjIHW < padj &
             abs(log2FCshrink_ashr) > lfc[1]) |>
    mutate(thresh_padj = padjIHW < padj & abs(log2FCshrink_ashr) > lfc[2])
  #
  df |>
    group_by(thresh_padjlfc == T) |> tally()
  df |>
    group_by(thresh_padj == T) |> tally()
  #
  df <- df |> dplyr::select(contains("thresh"), everything())
  n0 <- df |> filter(thresh_padjlfc == T) |> nrow()
  n00<- df |> filter(thresh_padj == T) |> nrow()
  df <- df |> mutate(genelabels = "")
  if (n00<=15){
    n1<-n00
    }else if (n0 > 15) {
    #n1 <- min(c(n0, 10))
    n1<-15
  }else if (n0 == 0 & n00!=0) {
    n1<-n00
  }
  
  df<-df |> arrange(padjIHW) 
    df$genelabels[1:n1] <- df$rowLabel[1:n1]
    df$genelabels[(n1+1):nrow(df)] <- ""
  
  
  xmax<-sort(abs(df$log2FCshrink_ashr),decreasing = T)[1]
  #volcano: https://hbctraining.github.io/DGE_workshop/lessons/06_DGE_visualizing_results.html
  df2<-df |> 
    select(log2FCshrink_ashr,thresh_padjlfc, thresh_padj, genelabels, padjIHW) |> 
    mutate(group = ifelse(thresh_padjlfc==T, paste0("p < 0.05;|lfc|>", lfc[1]), ifelse(thresh_padj==T, paste0("p < 0.05;|lfc|>", lfc[2]), ifelse(padjIHW>padj, paste0("p > ", padj), ""))))
  p1 <- ggplot(df, aes(x = log2FCshrink_ashr, y = -log10(padjIHW))) +
    xlim(-1*xmax, 1*xmax)+
    geom_point(
      data = df |> filter(thresh_padj == F),
      colour = okabe[9],
     # fill = okabe[9],
      size = 0.9,
      alpha = 0.3) +
    bb_theme()
  p1<-p1+ ylab(expression(-log10(P))) +
    xlab(expression(log2~(FC))) +
    ggtitle(cont)
   p1<- p1 + geom_point(data = df |> filter(thresh_padj == T),
      colour = okabe[2],
     # fill = okabe[2],
      size = 0.9,
      alpha = 0.5
    ) 
   p1<-p1+
    geom_point(
      data = df |> filter(thresh_padjlfc == T),
      colour = okabe[8],
      size = 0.9,
      alpha = 0.5
      ) 
    # geom_point(aes(colour = threshold_OE)) +
   
    p1<-p1+
      geom_text_repel(aes(label = genelabels),
                      col = 'black',
                      size = 1.5, seed=10, min.segment.length = 0,  max.overlaps = 15, segment.colour = "gray") 
   
    nL<-length(unique(df2$group))
    if(nL==2){
      cols<-c(okabe[2], okabe[9])
      }else{
        cols<-c(okabe[2], okabe[8], okabe[9])}
    p2<-df2 |> ggplot(aes(x = log2FCshrink_ashr, y = -log10(padjIHW), color = group)) + geom_point(aes(color=group), alpha=0.4) + 
      scale_color_manual(values = cols) + bb_theme() + theme(legend.position = "bottom")
    p2<-p2+xlim(-1*xmax, 1*xmax) 
    p2<- p2 + ylab(expression(-log10(P))) +
      xlab(expression(log2~(FC))) +
      ggtitle(cont) 
    p2<-p2 + geom_text_repel(aes(label = genelabels),
                         col = 'black',
                         size = 1.5,nudge_y = -0.001,max.overlaps = 15,min.segment.length = 0, seed=16,segment.colour = "gray") 
    p2<-p2 + geom_hline(yintercept = -log10(0.05), lty=2, linewidth = 0.2, colour = "gray")
    
    #geom_vline(xintercept = -1, lty=okabe[8]) +
    # scale_color_manual(c(okabe[9],okabe[8])) +
  
  #p1 <- p1 +
    #geom_vline(xintercept = -1*lfc,
    #           lty = 2,
    #           color = okabe[9]) +
    #geom_vline(xintercept = 1*lfc,
    #           lty = 2,
    #           color = okabe[9])
  
  #print(p1)
    print(p2)
  #ggsave(paste0(res_path, "figs/volcano_", cont, ext, ".pdf"), device =
   #        cairo_pdf)
  message("File saved to figs/volcano_", cont, ext, ".svg")
  ggsave(paste0(res_path, "figs/volcano_", cont, ext, ".svg"), width = 8)
  
}

#lapply
makeCont <- function(i = NULL) {
  comb_i <- combinations[, i]
  cond1 <- gs[comb_i[1]]
  cond2 <- gs[comb_i[2]]
  paste0(c("condition", cond2, cond1))
}

# PLOTS -------

## Heatmap Dist -----

MyPheatmapDist <- function(x,
                           fontsize = 12,
                           fontsize_row = 6,
                           fontsize_col = 6,
                           border_color = NA,
                           color = prgn,
                           cluster_cols = T,
                           cluster_rows = T,
                           clustering_distance_cols = NULL ,
                           clustering_distance_rows = NULL,
                           show_colnames = F,
                           main = "Add a title here") {
  pheatmap::pheatmap(
    x,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    border_color = border_color,
    color = color,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    clustering_distance_cols = clustering_distance_cols ,
    clustering_distance_rows = clustering_distance_rows,
    main = main,
    show_colnames = show_colnames
  )
}

## Heatmap ----
MyPheatmapClust <- function(x,
                            color = prgn,
                            kmeans_k = NA,
                            breaks = NA,
                            border_color = NA,
                            cellwidth = NA,
                            cellheight = NA,
                            scale = "row",
                            cluster_rows = TRUE,
                            cluster_cols = F,
                            clustering_distance_rows = "euclidean",
                            clustering_distance_cols = "euclidean",
                            clustering_method = "complete",
                            cutree_rows = NA,
                            cutree_cols = NA,
                            legend = TRUE,
                            legend_breaks = NA,
                            legend_labels = NA,
                            annotation_col = df2,
                            annotation_colors = ann_colors,
                            annotation_row = NA,
                            annotation = NA,
                            annotation_legend = TRUE,
                            annotation_names_row = TRUE,
                            annotation_names_col = TRUE,
                            drop_levels = TRUE,
                            show_rownames = F,
                            show_colnames = F,
                            main = "Add a title here",
                            fontsize = 12,
                            fontsize_row = 6,
                            fontsize_col = 5,
                            gaps_row = NULL,
                            gaps_col = NULL,
                            labels_row = NULL,
                            labels_col = NULL,
                            filename = NA,
                            width = NA,
                            height = NA,
                            silent = FALSE,
                            na_col = "#DDDDDD")
{
  pheatmap::pheatmap(
    x,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    border_color = border_color,
    color = color,
    scale = scale,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    main = main,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_colnames = show_colnames, 
    show_rownames = show_rownames,
    gaps_col = gaps_col,
    gaps_row = gaps_row,
    cutree_rows = cutree_rows,
    cutree_cols = cutree_cols,
    kmeans_k = kmeans_k
  )
}


## save heatmap ----
#from: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
save_pheatmap_png <- function(x,
                              filename,
                              width = 1200,
                              height = 1000,
                              res = 150) {
  png(filename,
      width = width,
      height = height,
      res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#find memory drags
# https://stackoverflow.com/questions/5749058/extend-memory-size-limit-in-r?rq=4
getsizes <- function() {
  z <- sapply(ls(envir = globalenv()), function(x)
    object.size(get(x)))
  mat<-as.matrix(as.matrix(rev(sort(z))))
  mat[,1]<-as.matrix(as.matrix(mat[,1]/(1000^2)))
  colnames(mat)<-"Mb"
  mat
 
  #(tmp <- as.matrix(rev(sort(z)))[1:10,]
}



#My Upset plot ----

MyUpsetPlot <- function(x = deseq_timept1,
                        mode = "distinct",
                        file = paste0(res_path, "figs/timePt1upsetDeseqDISTINCT.png"),
                        main = "TEST", 
                        png = T) {
  require(ComplexHeatmap)
  (m <- ComplexHeatmap::make_comb_mat(x, mode = mode))
  if(png == T){
  png(
    file,
    width = 14,
    height = 6,
    units = "in",
    res = 600
  )}
  if (mode == "intersect") {
    m <- m[comb_degree(m) > 1]
  }
  ss <- set_size(m)
  cs <- comb_size(m)
  ht <- UpSet(
    m,
    set_order = order(ss, decreasing = T),
    comb_order = order(comb_degree(m), -cs),
    #comb_col = RColorBrewer::brewer.pal(n = 11, "BrBG")[7:11][comb_degree(m)],
    #bg_col = "white",
    top_annotation = HeatmapAnnotation(
      "Intersection\nsize" = anno_barplot(
        cs,
        ylim = c(0, max(cs) * 1.1),
        border = FALSE,
        gp = gpar(fill = "black"),
        height = unit(4, "cm")
      ),
      annotation_name_side = "left",
      annotation_name_rot = 90
    ),
    left_annotation = rowAnnotation(
      "DEG set size" = anno_barplot(
        -ss,
        baseline = 0,
        axis_param = list(
          at = c(0, -1000, -2000),
          labels = c(0, 1000, 2000),
          labels_rot = 90,
          gp = gpar(fontsize = 8) # Add this line to reduce font size),
         # width = unit(3, "cm")
        ),
        border = FALSE,
        gp = gpar(fill = "black"),
        width = unit(4, "cm")
      ),
      set_name = anno_text(
        set_name(m),
        location = 0.5,
        just = "center",
        width = max_text_width(set_name(m)) + unit(2, "mm"),
        gp = gpar(fontsize = 9)
      )
    ),
    right_annotation = NULL,
    show_row_names = FALSE,
    use_raster = T
  )
  
  ht = draw(ht, column_title = paste0(main))
  od = column_order(ht)
  rod = row_order(ht)
  decorate_annotation("Intersection\nsize", {
    grid.text(
      cs[od],
      x = seq_along(cs),
      y = unit(cs[od], "native") + unit(2, "pt"),
      default.units = "native",
      just = c("left", "bottom"),
      gp = gpar(fontsize = 8, col = "#404040"),
      rot = 45
    )
  })
  
  decorate_annotation("DEG set size", {
    grid.text(
      ss[rod],
      x = unit(-ss[rod], "native") + unit(-3, "pt"),
      #
      #x = unit(-3, "native") + unit(20, "pt"),
      y = rev(seq_along(ss)),
      default.units = "native",
      just = c("centre", "bottom"),
      gp = gpar(fontsize = 8, col = "#404040"),
      rot = 90
    )
  })
  if(png == T){
  dev.off()
  }
}

#modified from gostplot ------
MyGoPlot <- function (gostres,
                      capped = F,
                      interactive = FALSE,
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
                      ),
                      thres = 0.05,
                      annot = T, 
                      title = "A title",
                      subtitle = "A subtitle")
{
  require(colorspace)
  if (is.null(pal)) {
    pal <- c(
      `GO:MF` = "#dc3912",
      `GO:BP` = "#ff9900",
      `GO:CC` = "#109618",
      KEGG = "#dd4477",
      REAC = "#3366cc",
      WP = "#0099c6",
      TF = "#5574a6",
      MIRNA = "#22aa99",
      HPA = "#6633cc",
      CORUM = "#66aa00",
      HP = "#990099"
    )
  }
  if (!("result" %in% names(gostres)))
    stop("Name 'result' not found from the input")
  if (!("meta" %in% names(gostres)))
    stop("Name 'meta' not found from the input")
  source_order <- logpval <- term_id <- opacity <- NULL
  term_size <- term_name <- p_value <- term_size_scaled <- NULL
  df <- gostres$result
  meta <- gostres$meta
  essential_names <- c("source_order",
                       "term_size",
                       "term_name",
                       "term_id",
                       "source",
                       "significant")
  if (!(all(essential_names %in% colnames(df))))
    stop(paste(
      "The following columns are missing from the result:",
      paste0(setdiff(essential_names, colnames(df)), collapse = ", ")
    ))
  if (!any(grepl("p_value", colnames(df))))
    stop("Column 'p_value(s)' is missing from the result")
  widthscale <- unlist(lapply(meta$query_metadata$sources, function(x)
    meta$result_metadata[[x]][["number_of_terms"]]))
  names(widthscale) <- meta$query_metadata$sources
  space <- 1000
  starts <- c()
  start <- 1
  starts[1] <- start
  if (!length(widthscale) < 2) {
    for (idx in 2:length(widthscale)) {
      starts[idx] <- starts[idx - 1] + space + widthscale[idx -
                                                            1]
    }
  }
  names(starts) <- names(widthscale)
  if (is.null(names(pal))) {
    names(pal) = meta$query_metadata$sources[1:length(pal)]
  }
  sourcediff = setdiff(meta$query_metadata$sources, names(pal))
  colors = grDevices::colors(distinct = TRUE)[grep(
    "gr(a|e)y|white|snow|khaki|lightyellow",
    grDevices::colors(distinct = TRUE),
    invert = TRUE
  )]
  if (length(sourcediff) > 0) {
    use_cols = sample(colors, length(sourcediff), replace = FALSE)
    pal[sourcediff] <- use_cols
  }
  if ("p_values" %in% colnames(df)) {
    p_values <- query <- significant <- NULL
    df$query <- list(names(meta$query_metadata$queries))
    df <- tidyr::unnest(data = df,
                        cols = c(p_values, query, significant))
    df <- dplyr::rename(df, p_value = p_values)
  }
  logScale <- function(input,
                       input_start = 1,
                       input_end = 50000,
                       output_start = 2,
                       output_end = 10) {
    m = (output_end - output_start) / (log(input_end) - log(input_start))
    b = -m * log(input_start) + output_start
    output = m * log(input) + b
    return(output)
  }
  xScale <- function(input,
                     input_start = 1,
                     input_end = sum(widthscale) +
                       (length(widthscale) - 1) * space,
                     output_start = 2,
                     output_end = 200) {
    m = (output_end - output_start) / (input_end - input_start)
    b = -m * input_start + output_start
    output = m * input + b
    return(output)
  }
    df$logpval <- -log10(df$p_value)
    df$opacity <- ifelse(df$significant, 0.9, ifelse(df$p_value == 
                                                       1, 0.1, 0.2))
    df$term_size_scaled = logScale(df$term_size)
    df <- df %>% dplyr::group_by(source) %>% 
      dplyr::mutate(order = xScale(source_order, 
input_start = 1, input_end = widthscale[source], output_start = starts[source], 
output_end = starts[source] + widthscale[source]))
    df$order <- xScale(df$order)
    if (capped) {
      df$logpval[df$logpval > 16] <- 17
      ymin <- -1
      ymax <- 18.5
      ticklabels <- c("0", "2", "4", "6", "8", "10", "12", 
                      "14", ">16")
      tickvals <- c(0, 2, 4, 6, 8, 10, 12, 14, 16)
    }else {
      ymin <- -1
      ymax <- ceiling(max(df$logpval)) + 5
      ticklabels <- ggplot2::waiver()
      tickvals <- ggplot2::waiver()
    }
    if (interactive) {
      sd <- crosstalk::SharedData$new(df, key = ~term_id)
    }else {
      sd <- df
    }
    sd <- sd |> 
      mutate(significant = ifelse(p_value <= thres, T, F))
    if(sum(grepl("highlighted", colnames(sd))) == 1){
    sd <- sd |> 
      mutate(Label = ifelse(highlighted == T, 
                            str_wrap(paste0(term_name, " (", term_id, ")"), width = 10), ""))
    sd <- sd |> 
      mutate(text = ifelse(highlighted == F, str_wrap(term_name, width = 10), ""))
    #if(sd |> filter(source == "KEGG") |> nrow() == 1){
      sd[sd$source == "KEGG",]$Label =  
        str_wrap(paste0(sd[sd$source == "KEGG",]$term_name, " (", 
        sd[sd$source == "KEGG",]$term_id, ")"), width = 10)
      sd[sd$source == "KEGG",]$text =  ""
    #}
}else{
  sd <- sd |> 
    mutate(Label = "", 
           text = ifelse(p_value <= thres, str_wrap(term_name, width = 10), ""))
}
    p <- ggplot2::ggplot(data = sd, 
                         ggplot2::aes(x = order, y = logpval, 
                         label = Label
                         ))
    p<-p+ggplot2::geom_point(data = sd[sd$significant ==F,], 
                             aes(color = "gray", 
                                 size = term_size_scaled, 
                                 alpha = opacity), 
                             show.legend = T)
    
    p<-p+ggplot2::geom_point(data=sd[sd$significant ==T,], 
                             ggplot2::aes(color = source, 
                                       size = term_size_scaled, 
                                       alpha = opacity), 
                               show.legend = TRUE) 
    p<-p+ 
      ggplot2::facet_wrap(~query, ncol = 1, scales = "free_x", 
                          shrink = FALSE) + ggplot2::ylab("-log10(p-adj)") 
      p<-p + ggtitle(label = title,subtitle = subtitle)
    
      p<-p+ 
        ggplot2::scale_color_manual(values = pal) 
      p<-p +
        ggplot2::theme(
          #legend.position = "none", 
        panel.border = ggplot2::element_blank(), 
        strip.text = ggplot2::element_text(size = 12, colour = "darkgrey"), 
strip.background = ggplot2::element_blank(),
axis.title.x = ggplot2::element_blank(), 
axis.text.x = ggplot2::element_text(
  #size = 8, 
angle = 45, hjust = 1), 
axis.ticks.x = ggplot2::element_blank(), 
axis.ticks.y = ggplot2::element_line(color = "grey", 
                                     linewidth = 0.5), 
axis.line.x = ggplot2::element_line(color = "grey", 
linewidth = 0), 
axis.line.y = ggplot2::element_line(linewidth = 0.4, 
color = "grey"))

#axis.title.y = ggplot2::element_text(size = 10, 
#margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))
p <- p + ylab(expression(-log[10]~p)) + xlab("")
p<-p+ 
  ggplot2::scale_alpha(range = c(0, 0.9), limits = c(0, 0.9)) + 
  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax), labels = ticklabels, breaks = tickvals) + 
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 
 205), breaks = (xScale(starts) + xScale(starts + widthscale))/2, labels = names(widthscale))
    
    for (s in names(widthscale)) {
      xstart = xScale(starts[s])
      xend = xScale(starts[s] + widthscale[s])
      p <- p + ggplot2::annotate("segment", x = xstart, xend = xend, 
                                 y = -1, yend = -1, size = 3, colour = pal[s])
    }
    
   
      #p<-p+ggplot2::geom_segment(aes(x = 0, xend = xend, y = -log10(thres), 
      #                           yend = -log10(thres)), size = 0.8, color = "grey")
      p<-p+ ggplot2::geom_hline(yintercept = -log10(thres),
                                linetype = "dashed", linewidth = 0.5, col = "gray")
   
p<-p+ 
  ggplot2::theme(plot.margin = ggplot2::margin(t = 2, r = 5, b = 20, l = 20, unit = "pt"))


    if(isTRUE(annot)){
     #p <- p + ggrepel::geom_text_repel(data = sd[sd$text!="",],aes(label = text), size = 2.5, nudge_x = -2, nudge_y = -2)
     p <- p + ggrepel::geom_label_repel(nudge_y = 4, nudge_x = 4, size = 2.5,)
    }
    if (capped) {
      #p <- p + ggplot2::annotate(geom = "text", x = 180, y = 16.2, 
       #                          label = "values above this threshold are capped", 
      #                           size = 2, color = "grey") + 
      #  ggplot2::geom_hline(yintercept = 16, 
      #linetype = "dashed", linewidth = 0.2, color = "grey")
    }

p <- p + bb_theme() 
    if (interactive) {
      p <- p |> plotly::ggplotly(tooltip = "text")
      p <- p |> plotly::highlight(on = "plotly_click", off = "plotly_doubleclick", 
                                   dynamic = FALSE, persistent = FALSE)
    }
    return(p)
  }
