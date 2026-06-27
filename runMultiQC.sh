#!/usr/bin/bash

base_dir="/home/bbitarello/scratch"
proj_dir="/home/bbitarello/projects/Ehux_Ehv201"

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
--title "After fastp run" \
--replace-names "${proj_dir}/data/sample-names-multiqc.tsv"
}

run_mqtc_fastp
          

# for salmon results:
out_dir=${base_dir}/multiqc/salmon-v1.10.3
mkdir -p ${out_dir}

cd ${base_dir}/salmon-v1.10.3

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

# for both
out_dir=${base_dir}/multiqc/fastp-and-salmon
mkdir -p ${out_dir}

cd ${base_dir}/

run_mqtc_fastp_and_salmon() {
multiqc  fastp-v024/ salmonQ-v1.10.3/ \
--force \
--outdir "${out_dir}" \
--zip-data-dir \
--export \
--verbose \
--title "After fastp and salmon quant runs" \
--replace-names "${proj_dir}/data/sample-names-multiqc.tsv"
}


#both
run_mqtc_fastp_and_salmon



