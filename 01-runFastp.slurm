#!/usr/bin/bash

#SBATCH -N 1-3
#SBATCH -c 4
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p high
#SBATCH --array=1-321
#SBATCH --job-name=FASTP
#SBATCH -o /home/bbitarello/scratch/logs/fastp-v024/fastp.%A-%a.log
#SBATCH -e /home/bbitarello/scratch/logs/fastp-v025/fastp.%A-%a.err

# PATHS
FASTQDIR="/home/bbitarello/scratch/FASTQ"
# OUTPUT
OUTDIR="/home/bbitarello/scratch/fastp-v024/"

mkdir -p $OUTDIR

# INPUT
METADATA="/home/bbitarello/projects/Ehux_Ehv201/data/metadata.tsv"
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $METADATA )
TIME=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $4 }' $METADATA )
RUNID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $5 }' $METADATA )
RUNLAB=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $6 }' $METADATA )
LANE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $7 }' $METADATA )
BARCODE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $8 }' $METADATA )





# Load environment and activate conda environment
source ~/.bashrc
eval "$(conda shell.bash hook)"
#conda activate fastpenv (fastp v 0.22)
#conda activate /home/bbitarello/miniconda3/envs/fastpenv
#conda activate fastpenv (fastp v 0.24)
conda activate /home/bbitarello/miniconda3/envs/fastp

# COMMAND

#fastp -v
fastp \
    -i ${FASTQDIR}/${RUNID}/${RUNID}_s_${LANE}_1_${BARCODE}.fastq.gz \
    -I ${FASTQDIR}/${RUNID}/${RUNID}_s_${LANE}_2_${BARCODE}.fastq.gz \
    -o ${OUTDIR}/${SAMPLE}_${RUNLAB}l${LANE}_r1.fp.fq.gz \
    -O ${OUTDIR}/${SAMPLE}_${RUNLAB}l${LANE}_r2.fp.fq.gz \
    -h ${OUTDIR}/${SAMPLE}_${RUNLAB}l${LANE}.html \
    -j ${OUTDIR}/${SAMPLE}_${RUNLAB}l${LANE}.json \
    --unpaired1 ${OUTDIR}/${SAMPLE}_${RUNLAB}l${LANE}_unpaired1 \
    --unpaired2 ${OUTDIR}/${SAMPLE}_${RUNLAB}l${LANE}_unpaired2 \
    --detect_adapter_for_pe


