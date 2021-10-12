eval "$(conda shell.bash hook)"
conda activate assembly
cd ~/data/tb
zcat sample1_1.fastq.gz | head -8
wc -l *
zcat sample1_1.fastq | head -2 | tail -1 | wc -c

cd ~/data/tb
spades.py -1 ~/data/tb/sample1_1.fastq.gz -2  ~/data/tb/sample1_2.fastq.gz -o sample1_asm
cd sample1_asm
quast -r ../tb.fasta -o quast contigs.fasta 
abacas.pl –h
cp ~/data/tb/tb.fasta .
abacas.pl -r tb.fasta -q contigs.fasta -p nucmer -b -d -a -m -N -g sample1 -o sample1_asm
bwa index -a is sample1_asm.fasta
bwa mem -k 20 -c 100 -L 20 -U 20 -M -T 50 sample1_asm.fasta ~/data/tb/sample1_1.fastq.gz ~/data/tb/sample1_2.fastq.gz | samtools sort -o sample1_asm.bam
samtools index sample1_asm.bam

cd ~/data/tb
spades.py -1 ~/data/tb/sample2_1.fastq.gz -2  ~/data/tb/sample2_2.fastq.gz -o sample2_asm
cd sample2_asm
quast -r ../tb.fasta -o quast contigs.fasta 
abacas.pl –h
cp ~/data/tb/tb.fasta .
abacas.pl -r tb.fasta -q contigs.fasta -p nucmer -b -d -a -m -N -g sample2 -o sample2_asm
bwa index -a is sample2_asm.fasta
bwa mem -k 20 -c 100 -L 20 -U 20 -M -T 50 sample2_asm.fasta ~/data/tb/sample2_1.fastq.gz ~/data/tb/sample2_2.fastq.gz | samtools sort -o sample2_asm.bam
samtools index sample2_asm.bam

cd ~/data/tb
spades.py -1 ~/data/tb/sample3_1.fastq.gz -2  ~/data/tb/sample3_2.fastq.gz -o sample3_asm
cd sample3_asm
quast -r ../tb.fasta -o quast contigs.fasta 
abacas.pl –h
cp ~/data/tb/tb.fasta .
abacas.pl -r tb.fasta -q contigs.fasta -p nucmer -b -d -a -m -N -g sample3 -o sample3_asm
bwa index -a is sample3_asm.fasta
bwa mem -k 20 -c 100 -L 20 -U 20 -M -T 50 sample3_asm.fasta ~/data/tb/sample3_1.fastq.gz ~/data/tb/sample3_2.fastq.gz | samtools sort -o sample3_asm.bam
samtools index sample3_asm.bam

cd ~/data/tb
samtools view -b sample1.bam Chromosome:79000-84000  | samtools fastq - > region.fastq

spades.py -s region.fastq -o region_assembly

