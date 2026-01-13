# Install salmon

Salmon: https://combine-lab.github.io/salmon/getting_started/#obtaining-salmon

conda config --add channels conda-forge

conda config --add channels bioconda

conda create -n salmon salmon

conda activate salmon


# Make index with decoys
./runDecoys.sh


# Make metadata
MakeMetadata.R

# Run salmon quantification
sbatch runSalmon.slurm

#sacct|grep 13668|grep COMPLETED|grep -v batch|wc

#multiqc: see multiqc dir
