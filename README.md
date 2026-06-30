# Ehux_Ehv201

Data processing was run on a Linux (Ubuntu) server with multiple cores. 

Some of the plots were ran in R on a Mac.

# 0: download genomic/annotation data and install everything
./0-get-downloads.sh

./0-install.sh


# fastp

##https://github.com/opengene/fastp

Rscript --vanilla 01-MakeMetadata.R

sbatch 01-runFastp.slurm

## ./check-slurm-status <jobid>

# Salmon

##https://github.com/COMBINE-lab/salmon

# Make index with decoys

./02-runDecoys.sh

# Run salmon quantification

sbatch 02-runSalmon.slurm



# 3. Multiqc

./03-runMultiQC.sh

# 4. Annotations

dge/04-annotations-host-new.Rmd
dge/04-annotations-virus-new.Rmd

# 4. DGE

##see dge/README.md





# Old

##KEGG pathway: https://rest.kegg.jp/link/ehx/pathway saved as KEGG-pathway.txt
