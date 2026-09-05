# Ehux_Ehv201

Data processing was run on a Linux (Ubuntu) server with multiple cores. 

Some of the plots were ran in R on a Mac.

# 0: download genomic/annotation data and install everything

```{bash}
./0-get-downloads.sh

./0-install.sh
```

# 1. Deal with metadata and run [fastp](https://github.com/opengene/fastp):

```{bash}
# make metadata for fastp runs

Rscript --vanilla 01-MakeMetadata.R

# map mislabelled samples

Rscript --vanilla 01-raw-to-processed-sample-names.R

# run fastp
sbatch 01-runFastp.slurm
```
Optional: If you're impatient too:

```{bash}
./check-slurm-status <jobid>
```

# 2. [Salmon](https://github.com/COMBINE-lab/salmon) analyses

Make index with decoys:

```{bash}
./02-runDecoys.sh
```

Run salmon quantification:

```{bash}
sbatch 02-runSalmon.slurm
```

# 3. [Multiqc](https://seqera.io/multiqc/) analyses

./03-runMultiQC.sh

# 4. Differential gene expression analyses

See dge/README.md


# 5. Plots for publication

After step 4 is done:


```{r}
Rscript --vanilla Plots-and-Tables-Publication.R
```



