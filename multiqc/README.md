
basepath="/home/bbitarello/scratch/"

# for fastp results:
cd ${basepath}/scratch
mkdir -p multiqc/fastp
conda activate py3.13
multiqc fastp --force -outdir multiqc/fastp

