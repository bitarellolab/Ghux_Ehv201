#!/opt/R/bin/Rscript
# Make sample subsets

# See also: 01-MakeMetadata.R (for salmon runs)

# Goals:

# 1. Remove weird symbols from sample names (like 1r, etc)
# 2. Fix scrambled salmon file names
# 3. Make ready to use salmon and metadata files.


## Load Packages ---

library(tidyverse)


## Paths ----
# base path
#(base_path <- path.expand("/Users/bbitarello/Library/CloudStorage/GoogleDrive-barbarabitarello@gmail.com/My Drive/rna-seq-host-virus/data_and_res_gh_repo/"))
# Salmon path
(base_path <- path.expand("~/Documents/GitHub/Ghux_Ehv201/"))
# adjust this to your actual file path
(base_path2 <- path.expand("~/Documents/Ehux-for-pub/salmonQ-v1.10.3/"))
#(base_path2 <- path.expand("~/Documents/Ehux-for-pub/salmonQ-v1.10.3"))


## Read in metadata (original) --------
sample_meta <- read_tsv(paste0(base_path,"data/sampleTable_80samp.txt"),
                        col_names = T,
                        col_types = cols())



## Fix weird names  ------------------
if(FALSE){ #only need to run once
sample_meta$names[which(grepl("r", sample_meta$names))] <-
  gsub("r", "", sample_meta$names[grepl("r", sample_meta$names)])

write_tsv(sample_meta, paste0(base_path,"data/sampleTable_80samp-fix.txt"))
}else{
sample_meta<-read_tsv(paste0(base_path,"data/sampleTable_80samp-fix.txt")) 
}

# 80 samp ------------------ -----------------
sampSet <- "80samp"

sample_meta80 <- sample_meta 
sample_meta80 <- tibble(sample_meta80)

write_tsv(sample_meta80, paste0(base_path,"data/sampleTable_", sampSet, ".txt"))

# Fix salmon files sample names
# Reason: some sequencing runs were labelled incorrectly. This script fixes that.


salmon_files80 <-
  file.path(paste0(base_path2,
                   sample_meta80$names,
                   "/quant.sf")) |>
  setNames(sample_meta80$names)
# salmon quant files labelled with these old names actually belong to the 
#samples listed in names_to_replace_new. E.g. sample HHQ_inf_t4_4  quant.sf is 
#saved as DMSO_cntl_t1_1/quant.sf

old_names <- names(salmon_files80)
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
names(salmon_files80) <- new_names

# check
length(salmon_files80) == nrow(sample_meta80)

# coldata for deseq
#reorder
salmon_files80 <- salmon_files80[sample_meta80$names]
coldata80 <- data.frame(salmon_files80, sample_meta80)
colnames(coldata80) <- gsub("salmon_files", "files", colnames(coldata80))

#check
names(salmon_files80)==sample_meta80$names

# save files
write_tsv(coldata80, paste0(base_path,"data/coldata_", sampSet, ".txt"))
#write_tsv(coldata, paste0(base_path,"data/coldata_", sampSet, "salmonv1.10.3.txt"))
saveRDS(salmon_files79, paste0(base_path, "data/salmon_files_", sampSet, ".rds"))
#saveRDS(salmon_files, paste0(base_path, "data/salmon_files_", sampSet, "salmonv1.10.3.rds"))

# For box folder

coldata80$files80<-gsub("/Users/bbitarello/Documents/Ehux-for-pub/", "Bitarello-FilesForPublication/raw-files/", coldata80$files80)
colnames(coldata80)[1]<-"files"
# save file
write_tsv(coldata80, paste0(base_path,"data/salmon-files-", sampSet, ".tsv"))



# 79samp  ------------------ ------------------
sampSet <- "79samp"
if(sampSet == "79samp"){
(samples_to_remove <- c(
  "HHQ_cntl_t1_3")) #see fastp results
}

sample_meta79 <- sample_meta |>
  dplyr::filter(!(names %in% samples_to_remove))

sample_meta79 <- tibble(sample_meta79)

write_tsv(sample_meta79, paste0(base_path,"data/sampleTable_", sampSet, ".txt"))

# Fix salmon files sample names
# Reason: some sequencing runs were labelled incorrectly. This script fixes that.


salmon_files79 <-
  file.path(paste0(base_path2,
            sample_meta79$names,
            "/quant.sf")) |>
  setNames(sample_meta79$names)
# salmon quant files labelled with these old names actually belong to the 
#samples listed in names_to_replace_new. E.g. sample HHQ_inf_t4_4  quant.sf is 
#saved as DMSO_cntl_t1_1/quant.sf

old_names <- names(salmon_files79)
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
names(salmon_files79) <- new_names

# check
length(salmon_files79) == nrow(sample_meta79)

# coldata for deseq
#reorder
salmon_files79 <- salmon_files79[sample_meta79$names]
coldata79 <- data.frame(salmon_files79, sample_meta79)
colnames(coldata79) <- gsub("salmon_files", "files", colnames(coldata79))

#check
names(salmon_files79)==sample_meta79$names

# save files
write_tsv(coldata79, paste0(base_path,"data/coldata_", sampSet, ".txt"))
#write_tsv(coldata, paste0(base_path,"data/coldata_", sampSet, "salmonv1.10.3.txt"))
saveRDS(salmon_files79, paste0(base_path, "data/salmon_files_", sampSet, ".rds"))
#saveRDS(salmon_files, paste0(base_path, "data/salmon_files_", sampSet, "salmonv1.10.3.rds"))


## 77samp  ------------------
sampSet<-"77samp"
if(sampSet == "77samp"){
  ## QC-d samples (77samp) -----
  
  (samples_to_remove <- c(
    "HHQ_cntl_t1_3",
    #see fastp results
    "DMSO_inf_t2_1",
    #see PCA with 79samp
    "DMSO_cntl_t2_1" #see PCA with 79samp
  )) 
}


sample_meta77 <- sample_meta |>
  dplyr::filter(!(names %in% samples_to_remove))

sample_meta77 <- tibble(sample_meta77)

write_tsv(sample_meta77, paste0(base_path,"data/sampleTable_", sampSet, ".txt"))

# Fix salmon files sample names
# Reason: some sequencing runs were labelled incorrectly. This script fixes that.


salmon_files77 <-
  file.path(paste0(base_path2,
            sample_meta77$names,
            "/quant.sf")) |>
  setNames(sample_meta77$names)
# salmon quant files labelled with these old names actually belong to the 
#samples listed in names_to_replace_new. E.g. sample HHQ_inf_t4_4  quant.sf is 
#saved as DMSO_cntl_t1_1/quant.sf

old_names <- names(salmon_files77)
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
names(salmon_files77) <- new_names

# check
length(salmon_files77) == nrow(sample_meta77)

# coldata for deseq
#reorder
salmon_files77 <- salmon_files77[sample_meta77$names]
coldata77 <- data.frame(salmon_files77, sample_meta77)
colnames(coldata77) <- gsub("salmon_files", "files", colnames(coldata77))

#check
names(salmon_files77)==sample_meta77$names

# save files
write_tsv(coldata77, paste0(base_path,"data/coldata_", sampSet, ".txt"))
#write_tsv(coldata, paste0(base_path,"data/coldata_", sampSet, "salmonv1.10.3.txt"))
saveRDS(salmon_files77, paste0(base_path, "data/salmon_files_", sampSet, ".rds"))
#saveRDS(salmon_files, paste0(base_path, "data/salmon_files_", sampSet, "salmonv1.10.3.rds"))

## 76samp -----

### Remove this extra one due to broken flask 

sampSet <- "76samp"

(samples_to_remove <- c(
  "HHQ_cntl_t1_3", #see fastp results
  "DMSO_inf_t2_1",#see PCA with 79samp. Outlier confirmed in two separate runs.
  "DMSO_cntl_t2_1", #see PCA with 79samp. Outlier confirmed in two separate runs.
  "DMSO_cntl_t1_1" # pers. comm. A. Platt//Remove this extra one due to broken flask 
))


sample_meta76 <- sample_meta76 |>
  dplyr::filter(!(names %in% samples_to_remove))
sample_meta76 |> filter(names %in% samples_to_remove)
sample_meta76 |> filter(names %in% samples_to_remove) #check

write_tsv(sample_meta76, paste0(base_path, "data/sampleTable_",sampSet, ".txt"))

# coldata for deseq
#reorder
salmon_files76 <- salmon_files77[sample_meta76$names]
names(salmon_files76)

coldata76 <- data.frame(salmon_files76, sample_meta76)
colnames(coldata76) <- gsub("salmon_files76", "files", colnames(coldata76))

# check
length(salmon_files76) == nrow(sample_meta76)
# coldata for deseq


#check
names(salmon_files76)==sample_meta76$names

# save files
write_tsv(coldata76, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files76, paste0(base_path, "data/salmon_files_", sampSet, ".rds"))
#saveRDS(salmon_files76, paste0(base_path, "data/salmon_files_", sampSet, "salmonv1.10.3.rds"))

## Part 2: Infected samples (39samp) -----

### set of samples

sampSet <- "39samp"

(inf_samples<- sample_meta76$names[grep("inf", sample_meta76$names)])

#fix meta
sample_meta39 <- sample_meta |>
    dplyr::filter(names %in% inf_samples)
write_tsv(sample_meta39, paste0(base_path, "data/sampleTable_",sampSet, ".txt"))

# coldata for deseq
#reorder
salmon_files39 <- salmon_files76[sample_meta39$names]
coldata39 <- data.frame(salmon_files39, sample_meta39)
colnames(coldata39) <- gsub("salmon_files39", "files", colnames(coldata39))

# check
length(salmon_files39) == nrow(sample_meta39)
# coldata for deseq
  

#check
names(salmon_files39)==sample_meta39$names

# save files
write_tsv(coldata39, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files39, paste0(base_path, "data/salmon_files_", sampSet, ".rds")) 

## Infected samples t2: t4 (29samp) -----

sampSet <- "29samp"


#fix meta
(sample_meta29 <-sample_meta39|>
  dplyr::filter(timePt != "t1"))

write_tsv(sample_meta29, paste0(base_path,"data/sampleTable_",sampSet, ".txt"))
# coldata for deseq
#reorder
salmon_files29 <- salmon_files39[sample_meta29$names]
coldata29 <- data.frame(salmon_files29, sample_meta29)
colnames(coldata29) <- gsub("salmon_files29", "files", colnames(coldata29))

# check
length(salmon_files29) == nrow(sample_meta29)
# coldata for deseq


#check
names(salmon_files29)==sample_meta29$names

# save files
write_tsv(coldata29, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files29, paste0(base_path, "data/salmon_files_", sampSet, ".rds")) 

## Infected samples t3: t4 (20samp) -----

sampSet <- "20samp"

#fix meta
(sample_meta20 <- sample_meta39|>
    dplyr::filter(timePt != "t1" & timePt !="t2"))

write_tsv(sample_meta20, paste0(base_path,"data/sampleTable_",sampSet, ".txt"))
# coldata for deseq
#reorder
salmon_files20 <-salmon_files39[sample_meta20$names]
coldata20 <- data.frame(salmon_files20, sample_meta20)
colnames(coldata20) <- gsub("salmon_files20", "files", colnames(coldata20))

# check
length(salmon_files20) == nrow(sample_meta20)
# coldata for deseq


#check
names(salmon_files20)==sample_meta20$names

# save files
write_tsv(coldata20, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files20, paste0(base_path, "data/salmon_files_", sampSet, ".rds")) 

## Infected samples (t1 only)(10samp) -----

### set of samples
sampSet <- "10samp"


#fix meta
(sample_meta10 <- sample_meta39|>
    dplyr::filter(timePt == "t4"))


# coldata for deseq
#reorder
salmon_files10 <- salmon_files39[sample_meta10$names]
coldata10 <- data.frame(salmon_files10, sample_meta10)
colnames(coldata10) <- gsub("salmon_files10", "files", colnames(coldata10))

# check
length(salmon_files10) == nrow(sample_meta10)
# coldata for deseq


#check
names(salmon_files10)==sample_meta10$names

# save files
write_tsv(coldata10, paste0(base_path, "data/coldata_", sampSet, ".txt"))
saveRDS(salmon_files10, paste0(base_path, "data/salmon_files_", sampSet, ".rds")) 

sessionInfo() %>%
  capture.output() %>%
  writeLines(paste0("logs/session-info-FixMetadata-",lubridate::today(),".txt"))

