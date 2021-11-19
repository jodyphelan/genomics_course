eval "$(conda shell.bash hook)"
conda install -n base -c conda-forge mamba -y
for file in ~/git/genomics_course/conda_env/*yaml; do 
	mamba env create -f $file
done
