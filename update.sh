# Any updates
#cd ~/data/tb
#mkdir bedaquiline
#cd bedaquiline
#wget https://tbdr.lshtm.ac.uk/static/files.tgz
#tar -xvf files.tgz
#rm files.tgz
##cd ~/Documents
##rm -rf ~/Documents/genomics_course/
##git clone https://github.com/jodyphelan/genomics_course.git
##cp ~/Documents/genomics_course/cheats/temp.krak ~/data/nanopore_activity/kraken/.temp.krak


#### Added as nanopore fix for 05/12/2022 course only. Future installs will not require this

source /home/user/miniconda3/etc/profile.d/conda.sh
source /home/user/miniconda3/etc/profile.d/mamba.sh
conda remove --name nanopore --all -y
wget https://raw.githubusercontent.com/jodyphelan/genomics_course/master/conda_env/nanopore.yaml
mamba env create -f nanopore.yaml -y
rm nanopore.yaml
sudo apt update
sudo apt install wget lsb-release
export PLATFORM=$(lsb_release -cs)
wget -O- https://cdn.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://cdn.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
sudo apt update
sudo apt install ont-guppy-cpu
conda activate nanopore
#mamba install -c bioconda minimap2 -y
#mamba install -c bioconda mafft -y
#mamba install -c bioconda bcftools -y
#mamba install -c bioconda samtools -y
#mamba install -c bioconda figtree -y 
cd /home/user/data/nanopore_activity/kraken/
retaxdump
