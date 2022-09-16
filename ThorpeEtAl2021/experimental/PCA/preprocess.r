## R version 4.1.1
## vcfR version 1.13.0

## extract the CDS from the core genome
## using the core genome in vcf format (obtained from snippy)
## and the annotation in gff format

library(vcfR)
vcf <- read.vcfR("core.vcf.gz")
gff <- read.table("26695.gff", sep="\t", quote="", nrows = 1713)

gffCDS <- gff[gff$V3 == "CDS",]
pos <- NULL
for(i in 1:nrow(gffCDS)) {
  pos <- c(pos, gffCDS$V4[i]:gffCDS$V5[i])
}
vcfb <- vcf[is.biallelic(vcf),]
vcfCDS <- vcfb[which(as.integer(vcfb@fix[,2]) %in% pos),]

write.vcf(vcfCDS, "coreCDS.vcf.gz")

