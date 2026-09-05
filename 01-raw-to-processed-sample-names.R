#!/opt/R/bin/Rscript
#################
# Goal of this script: make a file linking names of raw samples (FASTQ files) to correct sample names (e.g. DMSO_cntl_t3_1).
# Why this is needed: some samples were mislabelled by the sequencing center. Our scripts deal with that mislabelling when reading in the fastq files,
# but we also want to have a record for this for the purposes of sharing the raw data files with correct names.
# See also: 04-MakeMetadata.R (for salmon runs), 01-FixMetadata.R
# Note: analyses were performed reading in the mislabelled samples and correcting the mismatch within script.
#################
library(tidyverse)
library(readxl)
library(conflicted)
conflicts_prefer(dplyr::filter)

##
(base_path <- path.expand("~/Documents/GitHub/Ghux_Ehv201/"))
## adjust this to your actual file path
(base_path2 <- path.expand("~/Documents/Ehux-for-pub/salmonQ-v1.10.3/"))

## Read in metadata (original) 
meta_data<-read_tsv(paste0(base_path,"data/metadata.tsv")) 

names_to_replace_old <- c(
  "DMSO_cntl_t1_1",
  "DMSO_inf_t1_1",
  "HHQ_cntl_t2_2",
  "HHQ_inf_t2_2",
  "DMSO_cntl_t3_3",
  "DMSO_inf_t3_3",
  "HHQ_cntl_t4_4",
  "HHQ_inf_t4_4"
)

names_to_replace_new <- c(
  "HHQ_inf_t4_4",
  "HHQ_cntl_t4_4",
  "DMSO_inf_t3_3",
  "DMSO_cntl_t3_3",
  "HHQ_inf_t2_2",
  "HHQ_cntl_t2_2",
  "DMSO_inf_t1_1",
  "DMSO_cntl_t1_1"
)

df<-tibble(old=names_to_replace_old, new=names_to_replace_new)

names_to_replace_old <- c(
  "DMSO_cntl_t1_1",
  "DMSO_inf_t1_1",
  "HHQ_cntl_t2_2",
  "HHQ_inf_t2_2",
  "DMSO_cntl_t3_3",
  "DMSO_inf_t3_3",
  "HHQ_cntl_t4_4",
  "HHQ_inf_t4_4"
)

names_to_replace_new <- c(
  "HHQ_inf_t4_4",
  "HHQ_cntl_t4_4",
  "DMSO_inf_t3_3",
  "DMSO_cntl_t3_3",
  "HHQ_inf_t2_2",
  "HHQ_cntl_t2_2",
  "DMSO_inf_t1_1",
  "DMSO_cntl_t1_1"
)

df<-tibble(old=names_to_replace_old, new=names_to_replace_new) |> mutate(sample_name = old)
nrow(meta_data)
df2<-left_join(meta_data |> filter(sample_name %in% df$old), df) |> select(sample_name, old, new, everything())

df2<-df2 |> select(-c("sample_name", "old"))

df2<-df2 |> rename(sample_name = new)
#now the other part

df3<- meta_data |> filter(!(sample_name %in% df$old))|> select(sample_name,everything())
final<-bind_rows(df2, df3)
write_tsv(final, "data/metadata-corrected-labels.tsv")

