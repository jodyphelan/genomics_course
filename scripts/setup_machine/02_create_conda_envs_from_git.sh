eval "$(conda shell.bash hook)"
conda install -n base -c conda-forge mamba -y
mkdir ~/logs
conda create --name assembly --file ~/git/genomics_course/conda_env/assembly.explicit.txt 2> ~/logs/assembly.errlog > ~/logs/assembly.log
conda create --name eqtl --file ~/git/genomics_course/conda_env/eqtl.explicit.txt 2> ~/logs/eqtl.errlog > ~/logs/eqtl.log
conda create --name gwas --file ~/git/genomics_course/conda_env/gwas.explicit.txt 2> ~/logs/gwas.errlog > ~/logs/gwas.log
conda create --name mapping --file ~/git/genomics_course/conda_env/mapping.explicit.txt 2> ~/logs/mapping.errlog > ~/logs/mapping.log
conda create --name methylation --file ~/git/genomics_course/conda_env/methylation.explicit.txt 2> ~/logs/methylation.errlog > ~/logs/methylation.log
conda create --name microbiome --file ~/git/genomics_course/conda_env/microbiome.explicit.txt 2> ~/logs/microbiome.errlog > ~/logs/microbiome.log
conda create --name nanopore --file ~/git/genomics_course/conda_env/nanopore.explicit.txt 2> ~/logs/nanopore.errlog > ~/logs/nanopore.log
conda create --name phylogenetics --file ~/git/genomics_course/conda_env/phylogenetics.explicit.txt 2> ~/logs/phylogenetics.errlog > ~/logs/phylogenetics.log
conda create --name rnaseq --file ~/git/genomics_course/conda_env/rnaseq.explicit.txt 2> ~/logs/rnaseq.errlog > ~/logs/rnaseq.log
conda create --name variant_detection --file ~/git/genomics_course/conda_env/variant_detection.explicit.txt 2> ~/logs/variant_detection.errlog > ~/logs/variant_detection.log
conda create --name ml --file ~/git/genomics_course/conda_env/ml.explicit.txt 2> ~/logs/ml.errlog > ~/logs/ml.log