wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -b
~/miniconda3/condabin/conda init
source ~/.bashrc
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -y -c conda-forge mamba

cd ~/
git clone https://github.com/jodyphelan/genomics_course.git
cd ~/genomics_course//conda_env/
mamba env create -f assembly.yaml
mamba env create -f eqtl.yaml
mamba env create -f gwas.yaml
mamba env create -f mapping.yaml
mamba env create -f methylation.yaml
mamba env create -f microbiome.yaml
mamba env create -f nanopore.yaml
mamba env create -f phylogenetics.yaml
mamba env create -f rnaseq.yaml
mamba env create -f variant_detection.yaml

cd ~/bin/

ln -s ~/genomics_course/non_conda_programs/aliview/aliview .
ln -s ~/genomics_course/non_conda_programs/tablet/tablet .
ln -s ~/genomics_course/non_conda_programs/tempest/bin/tempest .

cd ~/
wget genomics.lshtm.ac.uk/data.tar.gz

tar -xvf data.tar.gz

