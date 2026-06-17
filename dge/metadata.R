#!/opt/R/bin/Rscript
# Make sample subsets

## See also: MakeMetadata.R (for salmon runs)

#Goals:

#1. Remove weird symbols from sample names (like 1r, etc)
#2. Fix scrambled salmon file names
#3. Make ready to use salmon and metadata files.


## Load Packages ---

library(tidyverse)
#base path
#(base_path <- path.expand("/Users/bbitarello/Library/CloudStorage/GoogleDrive-barbarabitarello@gmail.com/My Drive/rna-seq-host-virus/data_and_res_gh_repo/"))
# Salmon path
(base_path <- path.expand("~/Documents/GitHub/Ehux_Ehv201/"))
# adjust this to your actual file path
#(base_path2 <- path.expand("~/Documents/rosalind-rna-seq-mirror/"))
(base_path2 <- path.expand("~/Documents/Ehux-for-pub/salmonQ-v1.11.4/"))


## Read in metadata (original)
sample_meta <- read_tsv(paste0(base_path,"data/sampleTable_80samp.txt"),
                        col_names = T,
                        col_types = cols())


## Fix weird names  ------------------

sample_meta$names[which(grepl("r", sample_meta$names))] <-
  gsub("r", "", sample_meta$names[grepl("r", sample_meta$names)])

# Remove these samples from all analyses (based on QC) ------------------

(samples_to_remove <- c(
  "HHQ_cntl_t1_3",
  #see fastp results
  "DMSO_inf_t2_1",
  #see PCA with 79samp
  "DMSO_cntl_t2_1" #see PCA with 79samp
))



## QC-d samples (77samp) -----

sampSet <- "77samp"

#if(sampSet == "77samp"){
sample_meta <- sample_meta |>
  dplyr::filter(!(names %in% samples_to_remove))

sample_meta <- tibble(sample_meta)

write_tsv(sample_meta, paste0(base_path,"data/sampleTable_77samp.txt"))

# Fix salmon files sample names
# Reason: some sequencing runs were labelled incorrectly. This script fixes that.

salmon_files <-
  file.path(paste0(base_path2,
            sample_meta$names,
            "/quant.sf")) |>
  setNames(sample_meta$names)
# salmon quant files labelled with these old names actually belong to the 
#samples listed in names_to_replace_new. E.g. sample HHQ_inf_t4_4  quant.sf is 
#saved as DMSO_cntl_t1_1/quant.sf

old_names <- names(salmon_files)
new_names <- old_names
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

new_names[unlist(lapply(names_to_replace_old, function(x)
  which(old_names == x)))] <- names_to_replace_new

# assign
names(salmon_files) <- new_names

# check
length(salmon_files) == nrow(sample_meta)

# coldata for deseq
#reorder
salmon_files <- salmon_files[sample_meta$names]
coldata <- data.frame(salmon_files, sample_meta)
colnames(coldata) <- gsub("salmon_files", "files", colnames(coldata))

#check
names(salmon_files)==sample_meta$names

# save files

write_tsv(coldata, paste0(base_path,"data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files, paste0(base_path, "data/salmon_files_", sampSet, ".RDS"))

## QC-d samples (76samp) -----

### Remove this extra one due to broken flask 

(samples_to_remove2 <- c(
  "HHQ_cntl_t1_3",
  #see fastp results
  "DMSO_inf_t2_1",
  #see PCA with 79samp
  "DMSO_cntl_t2_1", #see PCA with 79samp
  "DMSO_cntl_t1_1" # pers. comm. A. Platt
))


sampSet <- "76samp"

#sample_meta77 <- read_tsv("data/sampleTable_77samp.txt", 
                         #show_col_types = FALSE)
salmon_files77 <- readRDS(paste0(base_path, "data/salmon_files_77samp.RDS"))

#fix meta
sample_meta76 <- read_tsv(paste0(base_path, "data/sampleTable_77samp.txt"), 
                          show_col_types = FALSE) |>
  dplyr::filter(!(names %in% samples_to_remove2))
write_tsv(sample_meta76, paste0("data/sampleTable_",sampSet, ".txt"))

# coldata for deseq
#reorder
salmon_files76 <- readRDS(paste0(base_path, "data/salmon_files_77samp.RDS"))[sample_meta76$names]
coldata76 <- data.frame(salmon_files76, sample_meta76)
colnames(coldata76) <- gsub("salmon_files76", "files", colnames(coldata76))

# check
length(salmon_files76) == nrow(sample_meta76)
# coldata for deseq


#check
names(salmon_files76)==sample_meta76$names

# save files
write_tsv(coldata76, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files76, paste0(base_path, "data/salmon_files_", sampSet, ".RDS"))

## Infected samples (39samp) -----

### set of samples

# Read in metadata (original)
sample_meta76 <- read_tsv(paste0(base_path, "data/sampleTable_76samp.txt"),
                        col_names = T,
                        col_types = cols())


(inf_samples<- sample_meta76$names[grep("inf", sample_meta76$names)])

sampSet <- "39samp"


#fix meta
sample_meta39 <- sample_meta76 |>
    dplyr::filter(names %in% inf_samples)
write_tsv(sample_meta39, paste0(base_path, "data/sampleTable_",sampSet, ".txt"))

# coldata for deseq
#reorder
salmon_files39 <- readRDS(paste0(base_path, "data/salmon_files_76samp.RDS"))[sample_meta39$names]
coldata39 <- data.frame(salmon_files39, sample_meta39)
colnames(coldata39) <- gsub("salmon_files39", "files", colnames(coldata39))

# check
length(salmon_files39) == nrow(sample_meta39)
# coldata for deseq
  

#check
names(salmon_files39)==sample_meta39$names

# save files
write_tsv(coldata39, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files39, paste0(base_path, "data/salmon_files_", sampSet, ".RDS")) 

## Infected samples t2: t4 (29samp) -----

sampSet <- "29samp"

#sample_meta39 <-read_tsv("data/sampleTable_39samp.txt")
#salmon_files39 <- readRDS("data/salmon_files_39samp.RDS")

#fix meta
(sample_meta29 <- read_tsv(paste0(base_path,"data/sampleTable_39samp.txt"))|>
  dplyr::filter(timePt != "t1"))

write_tsv(sample_meta29, paste0(base_path,"data/sampleTable_",sampSet, ".txt"))
# coldata for deseq
#reorder
salmon_files29 <- readRDS(paste0(base_path, "data/salmon_files_39samp.RDS"))[sample_meta29$names]
coldata29 <- data.frame(salmon_files29, sample_meta29)
colnames(coldata29) <- gsub("salmon_files29", "files", colnames(coldata29))

# check
length(salmon_files29) == nrow(sample_meta29)
# coldata for deseq


#check
names(salmon_files29)==sample_meta29$names

# save files
write_tsv(coldata29, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files29, paste0(base_path, "data/salmon_files_", sampSet, ".RDS")) 

## Infected samples t3: t4 (Xsamp) -----

sampSet <- "20samp"

#fix meta
(sample_meta20 <- read_tsv(paste0(base_path,"data/sampleTable_39samp.txt"))|>
    dplyr::filter(timePt != "t1" & timePt !="t2"))

write_tsv(sample_meta20, paste0(base_path,"data/sampleTable_",sampSet, ".txt"))
# coldata for deseq
#reorder
salmon_files20 <- readRDS(paste0(base_path, "data/salmon_files_39samp.RDS"))[sample_meta20$names]
coldata20 <- data.frame(salmon_files20, sample_meta20)
colnames(coldata20) <- gsub("salmon_files20", "files", colnames(coldata20))

# check
length(salmon_files20) == nrow(sample_meta20)
# coldata for deseq


#check
names(salmon_files29)==sample_meta29$names

# save files
write_tsv(coldata20, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files20, paste0(base_path, "data/salmon_files_", sampSet, ".RDS")) 

## Infected samples (t1 only)(10samp) -----

### set of samples
sampSet <- "10samp"


#fix meta
(sample_meta10 <- read_tsv(paste0(base_path, "data/sampleTable_39samp.txt")) |>
    dplyr::filter(timePt == "t4"))


# coldata for deseq
#reorder
salmon_files10 <- readRDS(paste0(base_path, "data/salmon_files_39samp.RDS"))[sample_meta10$names]
coldata10 <- data.frame(salmon_files10, sample_meta10)
colnames(coldata10) <- gsub("salmon_files10", "files", colnames(coldata10))

# check
length(salmon_files10) == nrow(sample_meta10)
# coldata for deseq


#check
names(salmon_files10)==sample_meta10$names

# save files
write_tsv(coldata10, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files10, paste0(base_path, "data/salmon_files_", sampSet, ".RDS")) 

sessionInfo() %>%
  capture.output() %>%
  writeLines(paste0("logs/sessionInfo_metadata ",date(),".txt"))

