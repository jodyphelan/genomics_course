eval "$(conda shell.bash hook)"
conda install -n base -c conda-forge mamba -y
ls ~/git/genomics_course/conda_env/*yaml | xargs -n1 -P4 \
	mamba env create -f
