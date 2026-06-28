library(data.table)
library(tidyverse)
source("dge/plot_funcs.R")
fastp <- fread(
  "~/Library/CloudStorage/OneDrive-brynmawr.edu/Ehux-project/rosalind-rna-seq-mirror/multiqc/fastp/multiqc_data_1/general_stats_table.tsv"
)
#fastp<-fastp |> pivot_longer(cols=colnames(fastp)[colnames(fastp)!="Sample"])

metadata <- fread("data/metadata.tsv")

fastp2 <- fastp |> separate(
  Sample,
  into = c("run_id", "lane", "samp", "read", "barcode"),
  sep = "_",
  remove = F
)
fastp2 <-
  fastp2 |> select(-lane)

fastp2 <- fastp2 |> rename(lane = samp)

fastp2 <- fastp2 |> select(-read)
fastp2 <- fastp2 |> mutate(lane = as.numeric(lane))

fastp3 <- tibble(left_join(fastp2, metadata))

fastp3 |> mutate(sample_name2 = paste0(run_id, "_l", lane, "_", sample_name)) |> select(Sample, sample_name2) |> write_tsv("data/sample-names-multiqc.tsv")
fastp4<-fastp3 |> mutate(sample_name2 = paste0(run_id, "_l", lane, "_", sample_name)) |>
  arrange(barcode)
  

fastp4|> select(Sample, sample_name2) |>
  write_tsv("data/sample-names-multiqc.tsv", col_names = F)
#fastp3 |> filter(sample_name=="HHQ_cntl_t1_3") |> select(-c(lane, barcode, run_id))
#"AGAGGCAACC-CTAATGATGG"

plot1<-fastp3 |> ggplot(aes(x = fct_reorder(sample_name, `Reads After Filtering`), y = `Reads After Filtering`)) + geom_col()  + ylab("Reads After Filtering") + xlab("Sample") + bb_theme() + theme(axis.text.y = element_text(size = 12), axis.text.x =element_text(size = 6, angle = 85,vjust = 1, hjust = 1.3))
plot1

# Salmon

read_tsv("~/")