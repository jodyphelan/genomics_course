# multiqc in microbiome env
eval "$(conda shell.bash hook)"
# conda activate microbiome
# pip install multiqc

# albacore in nanopore env
#conda activate nanopore
#pip install ~/git/genomics_course/non_conda_programs/albacore/*whl

# create symlinks to executables in ~/bin
mkdir ~/bin
echo "export PATH=\$PATH:~/bin" >> ~/.bashrc
# aliview
ln -s ~/git/genomics_course/non_conda_programs/aliview/aliview ~/bin/
# tempest
ln -s ~/git/genomics_course/non_conda_programs/tempest/bin/tempest ~/bin/
# tablet
ln -s ~/git/genomics_course/non_conda_programs/tablet/tablet ~/bin/
# artemis
ln -s ~/git/genomics_course/non_conda_programs/artemis/{art,act} ~/bin/
# tracer
ln -s ~/git/genomics_course/non_conda_programs/tracer/bin/tracer ~/bin/
# igv
ln -s ~/git/genomics_course/non_conda_programs/IGV_Linux_2.16.2/igv.sh ~/bin/igv


# mkdir ~/software/
# cd ~/software/
# wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.4.6_linux64.tar.gz
# tar -xvf ont-guppy-cpu_6.4.6_linux64.tar.gz
# rm ont-guppy-cpu_6.4.6_linux64.tar.gz
# echo "export PATH=\$PATH:~/software/ont-guppy-cpu/bin" >> ~/.bashrc
