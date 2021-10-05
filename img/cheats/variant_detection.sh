mkdir ~/data/tb/variants

cd ~/data/tb/variants/
bcftools mpileup -B -Q 23 -d 2000 -C 50 -f ~/data/tb/tb.fasta ~/data/tb/sample1.bam | bcftools call --ploidy 1 -mvO v - > sample1.raw.vcf
cat sample1.raw.vcf | vcfutils.pl varFilter -d 10 -D 2000 > sample1.filt.vcf

bcftools mpileup -B -Q 23 -d 2000 -C 50 -f ~/data/tb/tb.fasta ~/data/tb/sample2.bam | bcftools call --ploidy 1 -mvO v - > sample2.raw.vcf
cat sample2.raw.vcf | vcfutils.pl varFilter -d 10 -D 2000 > sample2.filt.vcf

gatk HaplotypeCaller -R ~/data/tb/tb.fasta -I ~/data/tb/sample1.bam -O sample1.gatk.raw.vcf -ploidy 1

bgzip sample1.filt.vcf
tabix -p vcf sample1.filt.vcf.gz
bgzip sample2.filt.vcf
tabix -p vcf sample2.filt.vcf.gz

cd ~/data/tb/variants/
delly call -o sample1.delly.bcf -q 20 -s 3 -g ~/data/tb/tb.fasta ~/data/tb/sample1.bam
delly call -o sample2.delly.bcf -q 20 -s 3 -g ~/data/tb/tb.fasta ~/data/tb/sample2.bam