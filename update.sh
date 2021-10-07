# Any updates
cd /home/user/Documents
git clone https://github.com/jodyphelan/genomics_course.git
mamba env create --force -f genomics_course/conda_env/assembly.yaml
mamba env create --force -f genomics_course/conda_env/nanopore.yaml

eval "$(conda shell.bash hook)"
conda activate nanopore
pip install /usr/local/src/ont_alba*

