# Intro
eval "$(conda shell.bash hook)"
conda activate day1
cd ~/data/malaria/
head -n 5 Pf3D7_05.fasta
zcat IT.Chr5_1.fastq.gz | head -n 12

# TB
cd ~/data/tb/
bwa index ~/data/tb/tb.fasta

## Sample 1
trimmomatic PE sample1_1.fastq.gz sample1_2.fastq.gz -baseout sample1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
bwa mem -c 100 -R "@RG\tID:sample1\tSM:sample1\tPL:Illumina" -M -T 50 ~/data/tb/tb.fasta sample1_1P sample1_2P | samtools view -b - | samtools sort -o sample1.bam -
samtools index sample1.bam
samtools flagstat sample1.bam
rm sample1_1P sample1_1U sample1_2P sample1_2U

## Sample 2
trimmomatic PE sample2_1.fastq.gz sample2_2.fastq.gz -baseout sample2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
bwa mem -c 100 -R "@RG\tID:sample2\tSM:sample2\tPL:Illumina" -M -T 50 ~/data/tb/tb.fasta sample2_1P sample2_2P | samtools view -b - | samtools sort -o sample2.bam -
samtools index sample2.bam
samtools flagstat sample2.bam
rm sample2_1P sample2_1U sample2_2P sample2_2U


# Malaria
cd ~/data/malaria/
bwa index ~/data/malaria/Pf3D7_05.fasta
bwa mem ~/data/malaria/Pf3D7_05.fasta ~/data/malaria/IT.Chr5_1.fastq.gz ~/data/malaria/IT.Chr5_2.fastq.gz | samtools view -b - | samtools sort -o IT.Chr5.bam -
samtools index IT.Chr5.bam
samtools view IT.Chr5.bam | grep "IL39_6014:8:61:7451:18170"


# Tradis
cd ~/data/tradis
zcat TraDIS_reads_chr1.fastq.gz | wc -l
perl remove_transposon_and_filter_FASTQ.pl ~/data/tradis/TraDIS_reads_chr1.fastq.gz 0  | gzip -c > TraDIS_reads_chr1.filt.fq.gz 
bowtie-build Bp.genome.fa Bp.genome.fa
bowtie -y -a -v 0 -S Bp.genome.fa -q TraDIS_reads_chr1.filt.fq.gz 2> TraDIS_reads_chr1.bowtieLog.txt | samtools view -b - | samtools sort -o TraDIS_reads_chr1.bam -
samtools index TraDIS_reads_chr1.bam
samtools view -b TraDIS_reads_chr1.bam chr1:90207-90668 | samtools flagstat -
