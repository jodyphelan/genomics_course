cd ~/data/methylation
ls
zcat tb_pb_14.ipd.csv | head
analyse_motif_ipd.py tb_pb_14.ipd.csv.gz CTCCAG tb_pb_14.assembly.fa
analyse_motif_ipd.py tb_pb_16.ipd.csv.gz CTCCAG tb_pb_16.assembly.fa 
plot_ipd.R tb_pb_14,tb_pb_16 CTCCAG CTCCAG.pdf
analyse_motif_ipd.py tb_pb_2.ipd.csv.gz GTAYNNNNATC tb_pb_2.assembly.fa
analyse_motif_ipd.py tb_pb_16.ipd.csv.gz GTAYNNNNATC tb_pb_16.assembly.fa
plot_ipd.R tb_pb_16,tb_pb_2 GTAYNNNNATC GTAYNNNNATC.pdf
ls *.motif_summary.csv > files.txt
combine_motifs.py files.txt unfiltered_motifs.csv
combine_motifs.py files.txt filtered_motifs.csv --min_qual 60


raxmlHPC -m GTRGAMMA -s pacbio.fasta -n pacbio.ML -p 11334 -k -f a -x 13243 -N 100


bcftools view pacbio.vcf.gz -H | wc -l
bcftools view pacbio.vcf.gz -H  -s tb_pb_10,tb_pb_14,tb_pb_11 -x | wc -l
bcftools view pacbio.vcf.gz -s tb_pb_10,tb_pb_14,tb_pb_11 -x | bcftools csq -f tb_genome.fasta -g tb_genome.gff  | grep Rv3263
bcftools view pacbio.vcf.gz -s tb_pb_10,tb_pb_14,tb_pb_11 -x | bcftools csq -f tb_genome.fasta -g tb_genome.gff  | bcftools query -f '[%POS\t%SAMPLE\t%TBCSQ{1}\n]' | grep Rv3263

