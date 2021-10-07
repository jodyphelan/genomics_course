eval "$(conda shell.bash hook)"
conda activate eqtl
cd ~/data/eqtl
R
library(MatrixEQTL)
library(ggplot2)

base.dir = "/home/user/data/eqtl/tables/"
snp_file_name = paste(base.dir, "snps_table.txt", sep="")
expression_file_name = paste(base.dir, "expression_table.txt", sep="")
output_file_name = "/home/user/data/eqtl/results_eqtl.txt"
output_file_name.cis = "/home/user/data/eqtl/results_eqtl_cis.txt"
expr = read.table(expression_file_name, sep="\t", header=T, row.names=1)
snp = read.table(snp_file_name, sep="\t", header=T, row.names=1)
head(expr)
head(snp)

pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 0
useModel = modelLINEAR

snpspos <- read.table(paste(base.dir, "snpspos.txt", sep=""), sep="\t", header=T)
genepos <- read.table(paste(base.dir, "genepos.txt", sep=""), sep="\t", header=T)

head(snpspos)
head(genepos)

errorCovariance = numeric()
snp = SlicedData$new()
snp$fileDelimiter = "\t"
snp$fileOmitCharacters = "NA"
snp$fileSkipRows = 1
snp$fileSkipColumns = 1
snp$fileSliceSize = 2000
snp$LoadFile(snp_file_name)
gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(expression_file_name)
cvrt = SlicedData$new()

me = Matrix_eQTL_main(
	snps = snp,
	gene = gene,
	cvrt = cvrt,
	output_file_name = output_file_name,
	pvOutputThreshold = pvOutputThreshold_tra,
	useModel = useModel,
	errorCovariance = errorCovariance,
	verbose = TRUE,
	output_file_name.cis = output_file_name.cis,
	pvOutputThreshold.cis = 1,
	snpspos = snpspos,
	genepos = genepos,
	cisDist = 200,
	pvalue.hist = TRUE,
	min.pv.by.genesnp = FALSE,
	noFDRsaveMemory = FALSE)

res <- me$cis$eqtls
head(res)

snps_table <- read.table("/home/user/data/eqtl/snps_eqtl_table.txt", sep="\t", header=T)
snps_table[snps_table$POS==2423785 | snps_table$POS==2424864,]

rv2162 <- read.table("/home/user/data/eqtl/rv2162c.txt", sep="\t", header=T)
rv2162

rv2162$lineage <- as.factor(rv2162$lineage)
bp <- ggplot(rv2162, aes(x=gene, y=expg, fill=lineage)) +
geom_boxplot() + geom_point(position=position_jitterdodge())
show(bp)

