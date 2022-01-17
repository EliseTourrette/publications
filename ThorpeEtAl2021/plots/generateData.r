## generate the RData object used to plot the graphs

###################################
## generate dNdSancestorSumSelCoeff
###################################

## change this variable 
pathwd <- "20201105_4_1"

d1 <- NULL
d2 <- NULL
for(i in 1:7) {
  d <- read.table(paste0("lastg1_",i,".txt"))
  d1 <- cbind(d1,d$V1[-c(1,2)])
  d <- read.table(paste0("lastg2_",i,".txt"))
  d2 <- cbind(d2,d$V1[-c(1,2)])
}
for(i in 8:14) {
  d <- read.table(paste0("lastg1_",i,".txt"))
  d1[,i-7] <- d1[,i-7]+d$V1[-c(1,2)]
  d <- read.table(paste0("lastg2_",i,".txt"))
  d2[,i-7] <- d2[,i-7]+d$V1[-c(1,2)]
}

data1 <- NULL
data2 <- NULL
selCoeff <- c(0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001)
for(i in 2:ncol(d1)) {
  data1 <- rbind(data1,data.frame(ratio = d1[,i]/d1[,1], neutral = d1[,1], s = selCoeff[i], pop = "pop1"))
  data2 <- rbind(data2,data.frame(ratio = d2[,i]/d2[,1], neutral = d2[,1], s = selCoeff[i], pop = "pop2"))
}
data <- rbind(data1,data2)
save(data, file = "dNdSancestor_20201105_4_1_g8000.RData")

## befor doing the for loop on i (ie the rec level)
## need to generate the first RData object by changing the pathwd variable in the first part of this section
for(i in c(1,3,4)) {
  pathwd <- paste0('20201105_',i,'_1')
  load(file = paste0('dNdSancestor_',pathwd,'_g8000.RData'))
  
  data$deleterious <- data$ratio * data$neutral
  
  data <- data.frame(deleterious = c(data$deleterious[data$s == -0.005] +
                                       data$deleterious[data$s == -0.002] +
                                       data$deleterious[data$s == -0.001] +
                                       data$deleterious[data$s == -0.0005] +
                                       data$deleterious[data$s == -0.0002] +
                                       data$deleterious[data$s == -0.0001]),
                     neutral = c(data$neutral[data$s == -0.005]),
                     pop = c(data$pop[data$s == -0.005]))
  data$ratio <- data$deleterious / data$neutral
  
  save(data, file = paste0("dNdSancestorSumSelCoeff_",pathwd,"_g8000.RData"))
}



###################################
## generate sampleData.RData
## and
## generate posMut_PATHWD.txt
## input for ancestryFullGenome.py
###################################

## !!!!!!!!!
## need to generate tag2_8000.txt by taking the last line of the file mutations_8000.txt
## not that the tag lines should also be removed from this file before using it on this part of the script (or do not forget to remove them before using the data)
# tag_8000.txt was generated from mutations_8000.txt -> 2nd part, after line with 'Tag:'
# tr -s ' '  '\n'< tag_8000.txt > tag2_8000.txt &


setwd("/Users/elise/Desktop/todaysWork")
library(ggplot2)

n0 <- 429936 ## twice the number of mutated sites (the tag file is separated in two: the first n0/2 lines represent the mutation id and the last n0/2 lines represent the generateion of migration from pop1 to pop2 for the mutation)
data <- read.table("tag2_8000.txt")
data0 <- data.frame(id = data$V1[1:(n0/2)], tag = data$V1[(n0/2+1):n0])
id0 <- data0$id[data0$tag == 8000]

## version of mutation file that is without the tags of the mutations that migrated from population 1 to population 2
data <- read.table("mutations_8000.txt", stringsAsFactor = FALSE, comment.char = "$")

data1 <- data[data$V4 == "p1",]
data2 <- data[data$V4 == "p2",]

data1 <- data1[order(data1$V7),]
data2 <- data2[order(data2$V7),]

n1 <- 10168 ## population size of population 1 (at generation 8000), can get it from the file out.csv
n2 <- 10059 ## population size of population 2 (at generation 8000)
data1$V12 <- data1$V12/n1
data2$V12 <- data2$V12/n2

data2 <- data2[-which(data2$V5 %in% id0 & data2$V12 == 1/n2),] ## remove the mutations that were migrated last in g = 8000 and with frequency 1/pop size

frc <- data2$V12[data2$V5 %in% data1$V5] - data1$V12[data1$V5 %in% data2$V5]
fr1 <- - data1$V12[!(data1$V5 %in% data2$V5)]
fr2 <- data2$V12[!(data2$V5 %in% data1$V5)]

fr <- c(frc,fr1,fr2)
type <- c(rep("in both pop", length(frc)), rep("only in non bottleneck pop", length(fr1)), rep("only in bottleneck pop", length(fr2)))
selc <- c(data2$V8[data2$V5 %in% data1$V5],data1$V8[!(data1$V5 %in% data2$V5)],data2$V8[!(data2$V5 %in% data1$V5)])
id <- c(data2$V5[data2$V5 %in% data1$V5],data1$V5[!(data1$V5 %in% data2$V5)],data2$V5[!(data2$V5 %in% data1$V5)])
pos <- c(data2$V7[data2$V5 %in% data1$V5],data1$V7[!(data1$V5 %in% data2$V5)],data2$V7[!(data2$V5 %in% data1$V5)])

sub <- data.frame(fr = fr, s = selc, type = type, id = id, pos = pos)
sub$mtype <- "neutral"
sub$mtype[sub$s < 0] <- "deleterious"

sub0 <- sub[sub$fr <= 0.15 & sub$fr >= -0.15,]
sub1 <- sub[sub$fr > 0.15 | sub$fr < -0.15,]
sub <- rbind(sub1, sub0[sample(1:nrow(sub0),1000),])

save(sub, file = "sampleData.RData")


write.table(unique(sub$pos), file = "posMut_20201105_4_2.txt", col.names = FALSE, row.names = FALSE)
