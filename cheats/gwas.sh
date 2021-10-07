eval "$(conda shell.bash hook)"
conda activate gwas
cd ~/data/gwas
plink --bfile MD --missing --out MD
plink --bfile MD --het --out MD
R CMD BATCH imiss-vs-het.Rscript
R CMD BATCH imiss_het_fail.Rscript
plink --bfile MD --indep-pairwise 200 5 0.5 --out MD
plink --bfile MD --extract MD.prune.in --genome --out MD
perl run-IBD-QC.pl MD
R CMD BATCH  plot-IBD.Rscript
R CMD BATCH plot-pca-results.Rscript
R CMD BATCH write_pca_fail.R
cat fail*txt | sort -k1 | uniq > fail_qc_inds.txt
plink --bfile MD --remove fail_qc_inds.txt --make-bed --out clean.MD
plink --bfile clean.MD --missing --out clean.MD
R CMD BATCH lmiss-hist.Rscript
plink --bfile clean.MD --test-missing --allow-no-sex --out clean.MD
perl run-diffmiss-qc.pl clean.MD
plink --bfile clean.MD --exclude fail-diffmiss-qc.txt --maf 0.01 --geno 0.05 --hwe 0.00001 --make-bed --out clean.final.MD
plink --bfile clean.final.MD --assoc --ci 0.95 --adjust --allow-no-sex --out final.MD.assoc
R CMD BATCH GWAS_plots.R
