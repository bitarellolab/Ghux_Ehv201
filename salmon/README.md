#install salmon

#https://combine-lab.github.io/salmon/getting_started/#obtaining-salmon
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n salmon salmon
conda activate salmon
