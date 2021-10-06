cd ~/data/tb
zcat sample1_1.fastq.gz | head -8
wc -l *
zcat sample1_1.fastq | head -2 | tail -1 | wc -c

##### Not sure!!
velveth sample1_k51 51 -fastq.gz -shortPaired ~/data/tb/sample1_1.fastq.gz ~/data/tb/sample1_2.fastq.gz
velvetg sample1_k51 -cov_cutoff 2 -ins_length 280 -ins_length_sd 80 -read_trkg no -min_contig_lgth 150 -exp_cov auto -scaffolding yes -min_pair_count 20 -unused_reads yes
VelvetOptimiser.pl --s 35 --e 51 -f '-fastq.gz -shortPaired ~/data/tb/sample1_1.fastq.gz ~/data/tb/sample1_2.fastq.gz'
#############

cd sample1_k45
cp ~/data/tb/tb.fasta .
abacas.pl -r tb.fasta -q contigs.fa -p nucmer -b -d -a -m -N -g sample1 -o k45


cd ~/data/tb/sample1_k45
bwa index -a is k45.fasta
bwa mem -k 20 -c 100 -L 20 -U 20 -M -T 50 k45.fasta ~/data/tb/sample1_1.fastq.gz ~/data/tb/sample1_2.fastq.gz | samtools sort -o k45.bam
samtools index k45.bam



cd ~/data/tb
samtools view -b sample1.bam Chromosome:79000-84000  | samtools fastq - > region.fastq
spades.py -s region.fastq -o region_assembly
