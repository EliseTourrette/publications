## R version 4.1.1

## NBR SITES

## for the GRAPES input, we need the number of syn and non syn SITES
## for that use the gff file as well as the globa alignment file
## 1) get the CDS from the gff file (the sequence from the reference)
## 2) use the global alignment to check whether the CDS are present in all strains or not. 
## If that is not the case, remove them (we are only interested in the core genome)
## 3) compute the number of synonymous and non synonymous sites
## based on Nei and Gojobori method
## ie for one codon, the nbr of syn sites = sum(1 to 3) fi, fi = proportion of syn changes for the ith nucleotide in the codon
## and the nbr of non syn sites is n = 3-s
## for the whole sequence, the nbr of syn sites = sum(codons) s
## and the nbr of non syn sites = 3C - S, C = nbr of codons


geneticCode = data.frame(codons = c( 'AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CAT','CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC','CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT','GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC','TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT','TTA','TTC','TTG','TTT'), aa = c('K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I','Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V','*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F'))

geneticCode$s <- 0
for(i in 1:nrow(geneticCode)) {
    codon0 <- unlist(strsplit(geneticCode$codons[i], split = ""))
    for(j in 1:3) {
        for(mut in c('A','T','C','G')) {
            if(codon0[j] != mut) {
                codon <- codon0
                codon[j] <- mut
                if(geneticCode$aa[i] == geneticCode$aa[geneticCode$codons == paste(codon, collapse = "")]) geneticCode$s[i] <- geneticCode$s[i] + 1
            }
        }
  }
}
geneticCode$s <- geneticCode$s/3
geneticCode$n <- 3-geneticCode$s
geneticCode$s[geneticCode$aa == '*'] <- 0
geneticCode$n[geneticCode$aa == '*'] <- 0

## load the sequences of the CDS 
## obtained via python
data <- unlist(read.table("sequence_CDS_ref.fa"))
data <- strsplit(data, split = "")
data <- lapply(data, function(x) {
                            y <- NULL; 
                            for(i in seq(1, length(x),3)) {
                                y <- c(y, paste(x[i:(i+2)], collapse = ""))
                            }
                            return(y);
                        }
                    )
data <- unlist(data)
data <- toupper(data)
td <- table(data)
tcodon <- attributes(td)$dimnames$data
## expected number of synonymous sites
LS <- sum(unlist(lapply(1:length(tcodon), function(x) td[x] * geneticCode$s[geneticCode$codons == tcodon[x]])))
## expected number of non synonymous sites
LN <- sum(unlist(lapply(1:length(tcodon), function(x) td[x] * geneticCode$n[geneticCode$codons == tcodon[x]])))


## NBR POLYMORPHISMS

library(PopGenome)
## 1) define the populations for all individuals
## 2) calculate the SFS
## 3) add the outgroup
## 4) calculate divergence with the outgroup

## load VCF file core.vcf, with gff file
data <- readData("vcf", format="VCF", gffpath="gff") 

pop <- read.table("pylori_populations.tsv", col.names = c("strain", "pop"))
pop$strain <- paste0("snippy_", pop$strain)
pop <- rbind(pop, c("Reference", "hspEuropeN")) ## if use vcf_outgroup, add the outgroup as a population
listPop <- lapply(unique(pop$pop), function(x) pop$strain[pop$pop == x])
names(listPop) <- unique(pop$pop)

data <- set.populations(data, listPop)

## set synonymous sites using a reference sequence
data <- set.synnonsyn(data, ref.chr="fasta/26695.fasta")

## SYNONYMOUS calculations
data.syn <- detail.stats(data, subsites = "syn")
data.syn <- neutrality.stats(data.syn, subsites = "syn")
data.syn <- diversity.stats(data.syn, subsites = "syn")
data.syn <- diversity.stats.between(data.syn, subsites = "syn")
MAF.syn <- data.syn@region.stats@minor.allele.freqs[[1]]
## if an outgroup has been set then returns the unfolded SFS
## if no outgroup has been set then returns the folded SFS
SFS.syn <- apply(MAF.syn, 1, table) 

## NON-SYNONYMOUS calculations
data.nonsyn <- detail.stats(data, subsites = "nonsyn")
data.nonsyn <- neutrality.stats(data.nonsyn, subsites = "nonsyn")
data.nonsyn <- diversity.stats(data.nonsyn, subsites = "nonsyn")
data.nonsyn <- diversity.stats.between(data.nonsyn, subsites = "nonsyn")
MAF.nonsyn <- data.nonsyn@region.stats@minor.allele.freqs[[1]]
SFS.nonsyn <- apply(MAF.nonsyn, 1, table) 



##########################################
## DIVERGENCE DATA
##########################################
## use the core VCF file 
## as well as the reference genome and the outgroup genome aligned on the ref
## for these two whole sequence, keep only the core CDS
## then find a way to caoncatenate the VCF and the "VCF outgroup-ref"
## that way, will be able to have the sites that are substitutions with Hp but polymorphic between Hp and outgroup

data <- readData("vcf_outgroup", format="VCF", gffpath="gff") 
dataSep <- data
## use separate vcf files for Hp core
## and outgroup vs Hp ref core
## then concatenate both vcf files
data <- concatenate.regions(data)

pop <- read.table("pylori_populations.tsv", col.names = c("strain", "pop"))
pop$strain <- paste0("snippy_", pop$strain)
pop <- rbind(pop, c("Reference", "hspEuropeN"), c("outgroup", "outgroup")) ## if use vcf_outgroup, add the outgroup as a population
listPop <- lapply(unique(pop$pop), function(x) pop$strain[pop$pop == x])
names(listPop) <- unique(pop$pop)

data <- set.populations(data, listPop)

## set synonymous sites using a reference sequence
data <- set.synnonsyn(data, ref.chr="fasta/26695.fasta")

data <- set.outgroup(data, c("outgroup")) ## only if use vcf_outgroup

dot <- MKT(data)
MKT <- as.data.frame(get.MKT(dot)[1,][[1]])
MKT$p1 <- sapply(strsplit(row.names(MKT), split = "/"), function(x) x[1])
MKT$p2 <- sapply(strsplit(row.names(MKT), split = "/"), function(x) x[2])

MKT <- MKT[MKT$p2 == "pop12",]

Pdiv_nonsyn <- MKT$P1_nonsyn
Pdiv_syn <- MKT$P1_syn

D_nonsyn <- MKT$D_nonsyn
D_syn <- MKT$D_syn



##########################################
## WRITE THE INPUT FILES
##########################################

## input files without divergence data
for(i in 1:(length(listPop)-1)) {
    write.table(names(listPop)[i], file = paste0(names(listPop)[i], ".folded.dofe"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(matrix(c("core_genes", length(data@populations[[i]]), LN, unname(SFS.nonsyn[[i]][-1]), LS, unname(SFS.syn[[i]][-1])),nrow = 1), file = paste0(names(listPop)[i], ".folded.dofe"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

## input files with divergence data
for(i in 1:(length(listPop)-1)) {
    write.table(names(listPop)[i], file = paste0(names(listPop)[i], ".folded.dofe"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(matrix(c("core_genes", length(dot@populations[[i]]), LN, unname(SFS.nonsyn[[i]][-1]), LS, unname(SFS.syn[[i]][-1]), LN, D_nonsyn[i], LS, D_syn[i]),nrow = 1), file = paste0(names(listPop)[i], ".folded.dofe"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}


