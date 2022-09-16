## R version 4.1.1

## create the VCF file for the divergent dataset (ie add the polymorphic sites to the Hp core VCF file)
vcf <- read.table("vcf/Hp_core.vcf", comment.char = '', skip = 3)
outgroup <- read.table("outgroup_core_diff.txt", skip = 1, sep = " ")
seq <- read.table("alignments/outgroupAlignment/Hacinonychis_aligned.aln", skip = 1)

outgroup <- outgroup[,-4]
colnames(outgroup) <- c("position", "ref", "outgroup")

## look at the sites that are present in Hp core VCF
seq <- strsplit(seq$V1, split = "")
seqVCF <- toupper(seq[[1]][as.numeric(vcf[-1, 2])])
vcfOUT <- cbind(vcf, c("outgroup",seqVCF))

colnames(vcfOUT) <- vcfOUT[1,]
vcfOUT <- vcfOUT[-1,]

## replace the nucleotide by the correct allele
## if the allele in the outgroup is the same as the ref -> 0
## if the allele is different -> 1
## plus, if the allele is different check whether it was already present in the alternative
## if that is not the case, add it, in upper case
for(i in 1:nrow(vcfOUT)) {
  if(vcfOUT$outgroup[i]  %in% vcfOUT$REF[i]) {
    vcfOUT$outgroup[i] <- "0"
  } else if((vcfOUT$outgroup[i] != ".") & (vcfOUT$outgroup[i] != "-") & !(vcfOUT$outgroup[i] %in% strsplit(vcfOUT$ALT[i], split = ",")[[1]])) {
        vcfOUT$ALT[i] <- paste(vcfOUT$ALT[i], toupper(vcfOUT$outgroup[i]), sep = ",")
        vcfOUT$outgroup[i] <- as.character(length(strsplit(vcfOUT$ALT[i], split = ",")[[1]]))
    } else if((vcfOUT$outgroup[i] != ".") & (vcfOUT$outgroup[i] != "-") & (vcfOUT$outgroup[i] %in% c(strsplit(vcfOUT$ALT[i], split = ",")[[1]]))) {
        vcfOUT$outgroup[i] <- which(strsplit(vcfOUT$ALT[i], split = ",")[[1]] == toupper(vcfOUT$outgroup[i]))
    }      
}


## add the sites that were not polymorphic in Hp
## but are when you compare the outgroup and Hp (in the core genome)
## first, remove the sites that are already present in VCFOUT
outgroup2 <- outgroup[!(outgroup$position %in% as.numeric(vcfOUT$POS)),]
outgroup2$ref <- toupper(outgroup2$ref)
outgroup2$outgroup <- toupper(outgroup2$outgroup)

## create the vcf
## assume that the sites that were not in the Hp core VCF are monomorphic for all Hp strains
vcfOUT2 <- as.data.frame(matrix(NA, ncol = ncol(vcfOUT), nrow = nrow(outgroup2)))
colnames(vcfOUT2) <- colnames(vcfOUT)
vcfOUT2$POS <- outgroup2$position
vcfOUT2$REF <- outgroup2$ref
vcfOUT2$ALT <- outgroup2$outgroup
vcfOUT2$outgroup <- "1"
vcfOUT2[,10:(ncol(vcfOUT2) - 1)] <- "0"
vcfOUT2[,1] <- unique(vcfOUT[,1])
vcfOUT2$ID <- unique(vcfOUT$ID)
vcfOUT2$QUAL <- unique(vcfOUT$QUAL)
vcfOUT2$FILTER <- unique(vcfOUT$FILTER)
vcfOUT2$INFO <- unique(vcfOUT$INFO)
vcfOUT2$FORMAT <- unique(vcfOUT$FORMAT)



## write vcfOUT as a standard vcf file
## need to drop " around the character strings
## and need to add the three lines that were previously skipped (will do it afterwards)
write.table(vcfOUT2, file = "HpOutgroup.vcf", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(vcfOUT, file = "coreOutgroup.vcf", quote = FALSE, row.names = FALSE, sep = "\t")

