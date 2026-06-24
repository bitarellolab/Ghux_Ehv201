# Packages ----
suppressPackageStartupMessages({
  library(extrafont)
  library(systemfonts)
  library(hrbrthemes)
  library(sysfonts)
  library(RColorBrewer)
  library(GGally)
  library(hrbrthemes)
  library(gcookbook)
  library(extrafont)
  library(showtext)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(ComplexHeatmap)
  library(gt)
  library(ggpubr)
  library(svglite)
  library(cols4all)
  library(circlize)
})
# Fonts ----
#extrafont::font_import(prompt = F)
#extrafont::loadfonts(quiet = T)
hrbrthemes::import_roboto_condensed()
sysfonts::font_add_google("Roboto")
sysfonts::font_add_google("Roboto Condensed")
gdtools::register_gfont(family = "Roboto Condensed")
systemfonts::get_from_google_fonts(family = "Roboto Condensed")
#fonts <- systemfonts::system_fonts()
systemfonts::fonts_as_import(family = "Roboto Condensed")
systemfonts::match_font("Roboto Condensed")
#extrafont::font_import(paths= "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/hrbrthemes/fonts/", recursive = T, prompt = F)

  # For tibble decimal points, no rounding
  old <- options(
    pillar.sigfig = 6,
    pillar.print_max = 5,
    pillar.print_min = 5,
    pillar.advice = FALSE
)

# Colors ----

plot_palette <- function(palette) {
  g <- ggplot2::ggplot(
    data = data.frame(
      x = seq_len(length(palette)),
      y = "1",
      fill = palette
    ),
    mapping = ggplot2::aes(
      x = x, y = y, fill = fill
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_void()
  return(g)
}

okabe <- c(
  "#000000",
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7",
  "#999999"
)

paired_cols <- c(
  "#A6CEE3",
  "#1F78B4",
  "#B2DF8A",
  "#33A02C",
  "#FB9A99",
  "#E31A1C",
  "#FDBF6F",
  "#FF7F00",
  "#CAB2D6",
  "#6A3D9A",
  "#FFFF99",
  "#B15928"
)



trt_vir_cols <- c(paired_cols[1:2], paired_cols[7:8])

#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#colors2 <-colorRampPalette( rev(brewer.pal(9, "Blues")) )(20)

#my_colors <- colorRampPalette(c("cyan", "deeppink3"))
htMapCols<-colorRampPalette(rev(viridisLite::magma(n=256,alpha = 1)))(100)[-c(90:100)]
#htMapCols <- 
l1<-length(htMapCols)
col_fun = colorRamp2(breaks = c(-3, 0, 3), colors = c(htMapCols[1],htMapCols[l1/2], htMapCols[l1]))
col_fun = colorRamp2(breaks = c(-3, 0, 3), colors = c(htMapCols[1],htMapCols[l1/2], htMapCols[l1]), space = "RGB")
### not run
#CyanPink <- colorRampPalette(c("cyan", "deeppink3"))
blues <- colorRampPalette(brewer.pal(9, "Blues"))(255)
blues2 <- colorRampPalette(brewer.pal(9, "Blues"))(20)
grays <- rev(hcl.colors(n = 5, palette = "grays"))
#hues <- rev(hcl.colors(10, "Magenta"))
#prgn <- colorRampPalette(rev(brewer.pal(name = 'PRGn', n = 11)[-c(6, 7)]))(255)
prgn <- colorRampPalette(rev(brewer.pal(name = 'PRGn', n = 10)))(255)
base::remove(l1); l1<-length(prgn);
col_fun_prgn<-colorRamp2(breaks = c(-3, 0, 3), colors = c(prgn[1],prgn[l1/2], prgn[l1]))
col_fun_prgn2<-colorRamp2(breaks = c(-3, 0, 3), colors = c(prgn[1],prgn[l1/2], prgn[l1]), space = "RGB")
prgn2 <- colorRampPalette(rev(brewer.pal(name = 'PRGn', n = 11)[-c(6, 7)]))(20)
prgn3<- colorRampPalette(rev(hcl.colors(palette = "Purple-Green", n = 11)))(255)

# Other ----
#showtext::showtext_opts(dpi = 150)

bb_theme <-  function() {
  theme(
    text = element_text(face = "plain", family = "Roboto Condensed"),
    axis.text.x = element_text(
      size = 13,
      family = "Roboto Condensed",
      face = "plain"
    ),
    axis.text.y = element_text(
      size = 13,
      family = "Roboto Condensed",
      face = "plain"
    ),
    axis.title.x = element_text(
      size = 14,
      hjust = 1,
      family = "Roboto Condensed",
      face = "plain"
    ),
    axis.title.y = element_text(
      size = 14,
      hjust = 1,
      family = "Roboto Condensed",
      face = "plain"
    ),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      family = "Roboto Condensed",
      face = "plain"
    ),
    plot.title = element_text(family = "Roboto Condensed", size = 16),
    strip.text =  element_text(size = 13),
    panel.background = element_blank(),
    panel.grid.major.x = element_line(color = "gray", linewidth = 0.1)
  )
}



#-----rstudio-font-settings
if(TRUE){
showtext::showtext_auto()
showtext::showtext_opts(dpi=300)

if (interactive()) {
  options(device = "RStudioGD")
} else {
  options(
    device = function(...)
      ragg::agg_png(...)
  )
}
}

