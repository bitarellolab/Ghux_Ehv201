library(tidyverse)
library(data.table)
source("dge/plot_funcs.R")
source("dge/Functions.R")
salmon <- fread("scratch-Aug2025/host/76samp/tables/salmon-output-76samp-host.tsv.gz")
annot <- fread(paste0("data/Annotations-host/annot-host-less-2026-07-21.tsv.gz"))
salmon <- salmon |>
  group_by(Name %in% annot$tx_id) |>
  rename(source = `Name %in% annot$tx_id`) |>
  ungroup() |>
  mutate(source = ifelse(source ==
    TRUE, "host", "virus"))


# reads per gene ----
# this is the filtered salmon file, containing only genes that survived filtering
salmon_filt <- fread("scratch-Aug2025/host/76samp/tables/salmon-output-processed-filt-76samp-host.tsv.gz")
# check if there are gene ids with more than one transcript in this set
length(unique(salmon_filt$ensembl_gene_id)) # 28,082
length(unique(salmon_filt$tx_id)) #28,090
# there are, so:
reads_per_gene_per_sample <- salmon_filt |>
  group_by(sample_name,ensembl_gene_id) |> 
  dplyr::summarise(across(c(n_reads, tpm), sum))
  
# reads per sample ----

reads_per_sample <- salmon |>
  group_by(sample_name) |>
  summarise(Reads_mapped = sum(NumReads))
summary((
  reads_per_sample |> 
    mutate(Reads_mapped = Reads_mapped / 1e6) |> 
    pull(Reads_mapped)
))

# reads mapped to host/virus per sample ----
reads_per_sample_per_source <- salmon |>
  group_by(sample_name, source) |>
  summarise(Reads_mapped = sum(NumReads))
(reads_per_sample_per_source2 <- reads_per_sample_per_source |>
  pivot_wider(names_from = source, values_from = Reads_mapped) |>
    mutate(sample_name = factor(sample_name)) |>
  mutate(Total = host + virus) |>
  mutate(p_host = host / Total, p_virus = virus / Total))

# sanity check ----

reads_per_sample_per_source2 |>
  filter(
    sample_name %in% c(
      "HHQ_cntl_t1_3",
      "DMSO_cntl_t1_1",
      "DMSO_inf_t2_1",
      "DMSO_cntl_t2_1"
    )
  )

# by time point ----
summary(reads_per_sample_per_source2 |> filter(grepl("inf", sample_name) & grepl("t1", sample_name)))
summary(reads_per_sample_per_source2 |> filter(grepl("inf", sample_name) & grepl("t2", sample_name)))
summary(reads_per_sample_per_source2 |> filter(grepl("inf", sample_name) & grepl("t3", sample_name)))
summary(reads_per_sample_per_source2 |> filter(grepl("inf", sample_name) & grepl("t4", sample_name)))


#normalized counts ----
normalizedC<-readxl::read_xlsx("~/Library/CloudStorage/Box-Box/RNA-Seq analysis resources/Bitarello-FilesForPublication/deseq_counts_contrastsOfInterest_76samp_host+annotation.xlsx", sheet = 5)
normalizedC<-normalizedC |> select(tx_id, locus_tag, DMSO_cntl_t1_2: HHQ_inf_t4_5)

normC_proc<-normalizedC |> 
  pivot_longer(cols=DMSO_cntl_t1_2: HHQ_inf_t4_5) |> 
  mutate(t = ifelse(grepl("t1", name), "t1", ifelse(grepl("t2", name), "t2", ifelse(grepl("t3", name), "t3", "t4")))) |>
  mutate(trt=factor(ifelse(grepl("DMSO", name), "DMSO", "HHQ")))


normC_proc |> 
  ggplot(aes(x = trt, y = value, group = trt)) +
  geom_violin() + geom_boxplot(width = 0.4)+facet_wrap(vars(t))
