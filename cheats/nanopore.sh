conda activate nanopore
cd ~/data/nanopore_activity/basecalling 
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5_reads -t 2 -s basecalling_output -o fastq --disable_filtering
cd ~/data/nanopore_activity/basecalling/basecalling_output/workspace
head fastq_runid_dcb79bbb002d7acdb262be3be0484eca7fc2f015_0.fastq
Rscript ~/data/nanopore_activity/scripts/MinIONQC.R -i ~/data/nanopore_activity/basecalling/basecalling_output/sequencing_summary.txt
cd ~/data/nanopore_activity/kraken
kraken --db ~/data/nanopore_activity/kraken/KDB/ --output temp.krak basecalled_reads.fastq
kraken-translate --db  ~/data/nanopore_activity/kraken/KDB/ temp.krak > kraken_output
cd ~/data/nanopore_activity/mapping
bwa index reference.fasta
bwa mem -x ont2d ./reference.fasta basecalled_reads.fastq > alignment.sam
samtools view -q 15 -b -S alignment.sam > alignment.bam
samtools sort alignment.bam -o sorted.bam
samtools index sorted.bam

cd ~/data/nanopore_activity/variant_calling 
samtools depth sorted.bam > depth_statistics

R
data<-read.table("depth_statistics")
plot(data$V3,type="l",xlab="Reference Position", ylab="read Depth")
quit()


bcftools mpileup -q 8 -B -I -Ou -f reference.fasta sorted.bam | bcftools call -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz
bcftools consensus -f reference.fasta calls.vcf.gz -o consensus_sequence.fasta


cd ~/data/nanopore_activity/phylogenetics
aliview zika_dataset.fasta
cat zika_dataset.fasta consensus_sequence.fasta > zika_all.fasta
mafft zika_all.fasta > zika_all_aligned.fasta
raxmlHPC -T 1 -m GTRGAMMA -s ./zika_all_aligned.fasta -n zika_phylogeny -p 11334 -k -f a -x 13243 -N 1000 -#1
cd ~/data/nanopore_activity/phylogenetics/figtree
figtree RAxML_bipartitions.zika_phylogeny