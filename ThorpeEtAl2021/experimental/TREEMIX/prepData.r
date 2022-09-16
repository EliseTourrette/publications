## R version 4.1.1

data <- read.table("core.noN.vcf.gz", comment.char = "", skip = 3)
sampleName <- unlist(data[1,-(1:10)])

data <- read.table("core.clust")
sampleName2 <- data$V1

j <- NULL
for(i in sampleName) j <- c(j, i %in% sampleName2)
  
j2 <- NULL
for(i in sampleName2) j2 <- c(j2, i %in% sampleName)

