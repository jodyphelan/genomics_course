# multiqc in microbiome env
eval "$(conda shell.bash hook)"
conda activate microbiome
pip install multiqc

# albacore in nanopore env
conda activate nanopore
pip install ~/git/genomics_course/non_conda_programs/albacore/*whl

# create symlinks to executables in ~/bin
mkdir ~/bin
# aliview
ln -s ~/git/genomics_course/non_conda_programs/aliview/aliview ~/bin/
# tempest
ln -s ~/git/genomics_course/non_conda_programs/tempest/bin/tempest ~/bin/
# tablet
ln -s ~/git/genomics_course/non_conda_programs/tablet/tablet ~/bin/
# artemis
ln -s ~/git/genomics_course/non_conda_programs/artemis/{art,act} ~/bin/
