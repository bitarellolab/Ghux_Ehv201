# ncbi datasets

##install ncbi_datasets-cli: https://github.com/ncbi/datasets
conda update -n base conda

conda create -n ncbi_datasets

conda activate ncbi_datasets

conda install -c conda-forge ncbi-datasets-cli


# fastp

conda config --add channels bioconda
conda create -n fastpenv
conda activate fastpenv
conda install -c bioconda fastp
fastp -v # 0.22.0


# multiqc

conda create --name py3.13 python=3.13
conda activate py3.13
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install bioconda::multiqc

multiqc #v1.35

## Salmon: 
#https://combine-lab.github.io/salmon/getting_started/#obtaining-salmon


conda activate py3.13

conda config --add channels conda-forge

conda config --add channels bioconda

conda install bioconda::salmon

salmon -v #1.11.4


