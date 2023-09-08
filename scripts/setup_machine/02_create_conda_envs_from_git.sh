eval "$(conda shell.bash hook)"
conda install -n base -c conda-forge mamba -y
conda create --name assembly --file ~/git/genomics_course/conda_env/assembly.explicit.txt
conda create --name eqtl --file ~/git/genomics_course/conda_env/eqtl.explicit.txt
conda create --name gwas --file ~/git/genomics_course/conda_env/gwas.explicit.txt
conda create --name mapping --file ~/git/genomics_course/conda_env/mapping.explicit.txt
conda create --name methylation --file ~/git/genomics_course/conda_env/methylation.explicit.txt
conda create --name microbiome --file ~/git/genomics_course/conda_env/microbiome.explicit.txt
conda create --name nanopore --file ~/git/genomics_course/conda_env/nanopore.explicit.txt
conda create --name phylogenetics --file ~/git/genomics_course/conda_env/phylogenetics.explicit.txt
conda create --name rnaseq --file ~/git/genomics_course/conda_env/rnaseq.explicit.txt
conda create --name variant_detection --file ~/git/genomics_course/conda_env/variant_detection.explicit.txt