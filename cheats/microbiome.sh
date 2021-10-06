cd fastq/
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" \
	> fastq_abs_paths
ls | grep -oE 'BB_[0-9]+' | sort -t _ -k 2 -n | uniq | \
	awk -v wd=$(pwd) 'OFS="\t" {print $1, wd "/" $1 "_1.fastq.gz", wd "/" $1 "_2.fastq.gz"}' \
	>> fastq_abs_paths
cd -
mv fastq/fastq_abs_paths .
time qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path fastq_abs_paths \
	--output-path fastq_imported.qza \
	--input-format PairedEndFastqManifestPhred33V2
time qiime demux summarize \
	--i-data fastq_imported.qza \
	--o-visualization fastq_imported.qzv
# get meta
cat meta.csv | tr ',' '\t' | cut -f 1,2 > meta.tsv

time qiime dada2 denoise-paired \
	--i-demultiplexed-seqs fastq_imported.qza \
	--p-trunc-len-f 190 \
	--p-trunc-len-r 190 \
	--p-n-threads 4 \
	--o-table table.qza \
	--o-representative-sequences rep_seqs.qza \
	--o-denoising-stats denoising_stats.qza

time qiime metadata tabulate \
	--m-input-file denoising_stats.qza \
	--o-visualization denoising_stats.qzv

time qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file meta.tsv 

time qiime feature-table tabulate-seqs \
  --i-data rep_seqs.qza \
  --o-visualization rep_seqs.qzv

rename 's/table/table_all_samples/' table*
time qiime feature-table filter-samples \
	--i-table table_all_samples.qza \
	--p-min-frequency 10000 \
	--o-filtered-table table.qza

time qiime feature-table summarize \
	--i-table table.qza \
	--o-visualization table.qzv \
	--m-sample-metadata-file meta.tsv 


time qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences rep_seqs.qza \
	--o-alignment aligned_rep_seqs.qza \
	--o-masked-alignment masked_aligned_rep_seqs.qza \
	--o-tree unrooted_tree.qza \
	--o-rooted-tree rooted_tree.qza

time qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted_tree.qza \
  --i-table table.qza \
  --p-sampling-depth 15000 \
  --p-n-jobs-or-threads 4 \
  --m-metadata-file meta.tsv \
  --output-dir core-metrics-results

time qiime feature-classifier classify-sklearn \
	--i-classifier db/gg-13-8-99-nb-classifier.qza \
	--i-reads rep_seqs.qza \
	--p-n-jobs 4 \
	--o-classification taxonomy.qza

time qiime metadata tabulate \
	--m-input-file taxonomy.qza \
	--o-visualization taxonomy.qzv


mkdir -p alpha_tests
time qiime diversity alpha-group-significance \
        --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
        --m-metadata-file meta.tsv \
        --o-visualization alpha_tests/faith_pd_group_significance.qzv

time qiime diversity alpha-group-significance \
        --i-alpha-diversity core-metrics-results/evenness_vector.qza \
        --m-metadata-file meta.tsv \
        --o-visualization alpha_tests/evenness_group_significance.qzv

time qiime diversity alpha-group-significance \
        --i-alpha-diversity core-metrics-results/shannon_vector.qza \
        --m-metadata-file meta.tsv \
        --o-visualization alpha_tests/shannon_group_significance.qzv

time qiime diversity alpha-group-significance \
        --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
        --m-metadata-file meta.tsv \
        --o-visualization alpha_tests/observed_features_group_significance.qzv


mkdir -p beta_tests

for metric in bray_curtis jaccard unweighted_unifrac weighted_unifrac; do
	qiime diversity beta-group-significance \
		--i-distance-matrix core-metrics-results/${metric}_distance_matrix.qza \
		--m-metadata-file meta.tsv \
		--m-metadata-column BV \
		--o-visualization beta_tests/${metric}_significance.qzv \
		--p-permutations 9999
done

mkdir -p rarefaction
time qiime diversity alpha-rarefaction \
	--i-table table.qza \
    --p-max-depth 27775 \
    --m-metadata-file meta.tsv \
    --o-visualization rarefaction/alpha_rarefaction.qzv


time qiime taxa barplot \
	--i-table table.qza \
	--i-taxonomy taxonomy.qza \
	--m-metadata-file meta.tsv \
	--o-visualization taxa_barplot.qzv

