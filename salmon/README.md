# Make index with decoys
./01-runDecoys.sh


# Make metadata
MakeMetadata.R

# Run salmon quantification
sbatch runSalmon.slurm

#sacct|grep 13668|grep COMPLETED|grep -v batch|wc

#multiqc: see multiqc dir
