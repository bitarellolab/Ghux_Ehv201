# Ehux_Ehv201

# 0: download genomic/annotation data
bash 0-get-downloads.sh

##also: #KEGG pathway: https://rest.kegg.jp/link/ehx/pathway

# Salmon

salmon/README.md

## Install salmon

##Salmon: https://combine-lab.github.io/salmon/getting_started/#obtaining-salmon

conda config --add channels conda-forge

conda config --add channels bioconda

conda create -n salmon salmon

conda activate salmon

## prep salmon decoys

bash 01-runDecoy.sh

## run salmon


# 3. DGE



