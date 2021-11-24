# Any updates
conda activate base
wget https://raw.githubusercontent.com/jodyphelan/genomics_course/master/conda_env/nanopore.yaml
mamba env create -f nanopore.yaml -y
rm nanopore.yaml
