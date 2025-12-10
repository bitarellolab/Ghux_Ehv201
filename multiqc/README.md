#install multiqc
conda create --name py3.13 python=3.13
conda activate py3.13
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install multiqc
