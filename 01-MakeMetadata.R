#!/opt/R/bin/Rscript
library(tidyverse)
library(readxl)
library(conflicted)
conflicts_prefer(dplyr::filter)

# Goals: make metadata tsv file. See also: FixMetadata.R

study_info <- 
  read_excel("data/AAA-StudyInfo.xls") |>
  filter(RULA_Run %in% c("FGC2634", "FGC2639"))

#make metadata
meta_data <- 
  study_info |>
  dplyr::select(COND_name, SAMP_sid, RULA_Run, RULA_Lane, RULA_Barcode) |>
  separate(COND_name, c("treatment", "virus", "timept"), sep = " ") |>
  mutate(timept = tolower(timept),
         timept = str_replace(timept, "tp", "t"),
         virus = recode(virus, "Control" = "cntl", "Infected" = "inf"),
         SAMP_sid = str_extract(SAMP_sid, "(\\d)$"),
         run_label = sub("^(F)GC\\d{2}(\\d{2})$", "\\1\\2", RULA_Run)) |>
  unite("condition", c(treatment, virus, timept), sep = "_", remove = FALSE) |>
  unite(sample_name, c(condition, SAMP_sid), sep = "_", remove = FALSE) |>
  select(sample_name, treatment, virus, timept, 
         run_id = RULA_Run, run_label, lane = RULA_Lane, barcode = RULA_Barcode)

meta_data


write_tsv(meta_data, "data/metadata.tsv")


sessionInfo() %>%
  capture.output() %>%
  writeLines(paste0("logs/session-info-FixMetadata-",lubridate::today(),".txt"))
