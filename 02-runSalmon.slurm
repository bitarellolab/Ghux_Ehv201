#!/bin/bash
#SBATCH --job-name=salmon-quant
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=16
#SBATCH -p high
#SBATCH --array=1-80
#SBATCH -o /home/bbitarello/scratch/final/logs/2026/salmon-quant.%A-%a.log
#SBATCH -e /home/bbitarello/scratch/final/logs/2026/salmon-quant.%A-%a.err


# INPUT and OUTPUT
DECOYPATH="/home/bbitarello/scratch/salmonI/decoyIndex"
# OUTPUT (change to decoy2 when using that)
OUTDIR="/home/bbitarello/scratch/salmon-quant-2026"
FILEEND=".fp.fq.gz"

mkdir -p ${OUTDIR}
FFASTQDIR="/home/bbitarello/scratch/fastp"
METADATA="/home/bbitarello/scratch/sample_ids.tsv"
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $METADATA )
cat $SAMPLE
#Load environment and activate conda environment
source ~/.bashrc
eval "$(conda shell.bash hook)"
#conda activate salmon
conda activate /home/bbitarello/miniconda3/envs/py3.13


run_salmon() {
salmon quant \
--threads 20 \
--index "${DECOYPATH}" \
--libType A \
--validateMappings \
-1 \
${FFASTQDIR}/${SAMPLE}_F34l1_r1${FILEEND} \
${FFASTQDIR}/${SAMPLE}_F34l2_r1${FILEEND} \
${FFASTQDIR}/${SAMPLE}_F39l1_r1${FILEEND} \
${FFASTQDIR}/${SAMPLE}_F39l2_r1${FILEEND} \
-2 \
${FFASTQDIR}/${SAMPLE}_F34l1_r2${FILEEND} \
${FFASTQDIR}/${SAMPLE}_F34l2_r2${FILEEND} \
${FFASTQDIR}/${SAMPLE}_F39l1_r2${FILEEND} \
${FFASTQDIR}/${SAMPLE}_F39l2_r2${FILEEND} \
-o \
${OUTDIR}/${SAMPLE}
}

run_salmon

