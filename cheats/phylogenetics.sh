eval "$(conda shell.bash hook)"
conda activate phylogenetics
cd ~/data/phylogenetics/
raxmlHPC -m GTRGAMMA -s H1N1.flu.2009.fas -n H1N1.flu.2009.ML -p 11334 -k -f a -x 13243 -N 100
mv RAxML_bipartitions.H1N1.flu.2009.ML RAxML_bipartitions.H1N1.flu.2009.ML.tre
