#install fastp: https://github.com/OpenGene/fastp

conda config --add channels bioconda
conda create -n fastpenv
conda activate fastpenv
conda install -c bioconda fastpenv


#run
sbatch run_fastp.slurm

#sacct|grep 13347|grep COMPLETED|grep -v batch|wc


#multiqc: see multiqc dir
