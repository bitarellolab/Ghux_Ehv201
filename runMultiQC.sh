#!/usr/bin/bash

base_dir="/home/bbitarello/scratch"

#Load environment and activate conda environment
source ~/.bashrc
eval "$(conda shell.bash hook)"
#conda activate salmon
conda activate /home/bbitarello/miniconda3/envs/py3.13

# for fastp results:
out_dir=${base_dir}/multiqc/fastp-v024
mkdir -p ${out_dir}

cd ${base_dir}/fastp-v024

run_mqtc_fastp() {
multiqc . \
--force \
--outdir "${out_dir}" \
--zip-data-dir \
--export \
--verbose \
--title "After fastp run"
}

#run_mqtc_fastp
          

# for salmon results:
out_dir=${base_dir}/multiqc/salmonQ-v1.11.4
mkdir -p ${out_dir}

cd ${base_dir}/salmonQ-v1.11.4

run_mqtc_salmon() {
multiqc . \
--force \
--outdir "${out_dir}" \
--zip-data-dir \
--export \
--verbose \
--title "After salmon quant run"
}

run_mqtc_salmon


