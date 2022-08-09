## this script will analyse and plot the ancestry matrix 
## obtained by the python script ancestry.py (or ancestry2.py, 
## depending if want to look at only the polymorphic sites or at all sites)
## (because I am useless at handling matrix and plotting in python...)

# setwd("/home/elise/Desktop/WORK/articleHP/simu10/simu10/analysesAncestry")
library(ggplot2)
library(tidyr)


## function to calculate the number of blocks
## via the size (in bp) of each block
## for one individual
## also return the type of the block
nbBlock <- function(data){
    nbBlock <- 1
    type <- data[1]
    begin <- 1
    lengthBlock <- NULL
    for(i in 2:length(data)) {
        if(data[i] != data[i-1]) {
            nbBlock = nbBlock +1
            type <- c(type, data[i])
            lengthBlock <- c(lengthBlock, i-begin)
            begin <- i
        }
    }
    lengthBlock <- c(lengthBlock, length(data)-begin)
    return(matrix(c(type, lengthBlock), byrow = TRUE, nrow = 2))
}

## function that to the correspondance between
## the number used in the python analysis
## in their correspondence in character string
renameAncestry <- function(data) {
    data$ancestry[data$ancestry == 1] <- "\nnon-bottleneck \npre admixture \n"
    data$ancestry[data$ancestry == 2] <- "\nbottleneck \npre admixture \n"
    data$ancestry[data$ancestry == 3] <- "\nnon-bottleneck \npost admixture \n"
    data$ancestry[data$ancestry == 4] <- "\nbottleneck \npost admixture \n"
    data$ancestry[data$ancestry == 5] <- "\nbefore split \n"
    data$ancestry[data$ancestry == 0] <- "\nunknown \n"
    return(data)
}

rec <- "noRec"
G <- 8000
interval <- 100

ancestryString <- c("\nnon-bottleneck \npre admixture \n", "\nbottleneck \npre admixture \n", "\nnon-bottleneck \npost admixture \n", "\nbottleneck \npost admixture \n", "\nbefore split \n", "\nunknown \n")
ancestryColor <- c("red", "blue", "orange", "purple", "green4", "white")

ancestryList <- NULL
ixA <- 1

for(file in system("ls ancestryMatrix", intern = TRUE)) {

rec <- strsplit(file, split = "_")[[1]][1]
G <- as.numeric(strsplit(file, split = "_")[[1]][3])

## load the ancestry matrix
## format pathwd_ancestryMatrix_G_cutoff.txt
data <- read.table(paste0("./ancestryMatrix/", rec,"_ancestryMatrix_", G,"_1000.txt"))

## load the vector of positions
## format pathwd_positions_G_cutoff.txt
posMut <- read.table(paste0("./positions/", rec,"_positions_", G,"_1000.txt"))
posMut <- unlist(posMut)[-c(1:3)]

## make sure that the data is numeric
## to be able to plot it afterwards
data2 <- apply(data,1,as.numeric)

## some data checks
#unlist(apply(data, 1, table))
#heatmap(data2, Rowv = NA)

## transpose the matrix
## to get row = individual and column = positions
data2 <- t(data2)

block <- apply(data2, 1, nbBlock)
block2 <- matrix(unlist(block), ncol = 2, byrow = TRUE)
block2[,2] <- block2[,2]*interval
block2 <- data.frame(ancestry = block2[,1], size = block2[,2])
## Attention: the number of locus for a block is not necessarily equal to its size
## it will depend on the size of he interval chosen
## for example if look at every site
## or if look at a site every 100 bp
## can get the interval size by looking at the position vector

## calculate the average size of the blocks
## depending on their ancestry
tmp <- aggregate(block2["size"], by = list(ancestry = block2$ancestry), mean)
meanSize <- data.frame(ancestry = c(1:5,0), size = rep(NA,6))
for(i in meanSize[,1]) {meanSize[ifelse(i == 0,6,i),2] <- ifelse(i %in% tmp[,1], tmp[tmp[,1] == i,2], NA)}
meanSize <- renameAncestry(meanSize)

data2 <- data.frame(ind = 1:nrow(data2), data2)
colnames(data2) <- c("ind", posMut)
indDendro <- as.dendrogram(hclust(d = dist(x = data2[,-1])))
indOrder <- order.dendrogram(indDendro)

data3 <- data2 %>% pivot_longer(!ind, names_to = "position", values_to = "ancestry")
data3 <- renameAncestry(data3)
data3$ind <- factor(x = data3$ind, levels = data2$ind[indOrder], ordered = TRUE)
data3$ancestry <- factor(x = data3$ancestry, levels = ancestryString, ordered = TRUE) 

ggplot(data3) +
  geom_tile(aes(x = position, y = ind, fill = as.factor(ancestry))) +
  scale_fill_manual(name = "ancestry \nmean block size", breaks = ancestryString, labels = paste0(ancestryString, round(meanSize$size), "\n"), values = ancestryColor)
ggsave(paste0("./plots/", rec, "_", G, ".jpeg"), width = 16, height = 9)


ancestryList[[ixA]] <- list(rec, G, meanSize)
ixA <- ixA + 1

    
}


save(ancestryList, file = "ancestry.RData")

ancestryData <- data.frame(rec = NULL, G = NULL, ancestry = NULL, size = NULL)

for(i in 1:length(ancestryList)) {
    anc <- ancestryList[[i]]
    ancestryData <- rbind(ancestryData, data.frame(rec = anc[[1]], G = anc[[2]], ancestry = anc[[3]]$ancestry, size = anc[[3]]$size))
}



ggplot(ancestryData) +
  geom_point(aes(x = G, y = size, color = ancestry)) +
  scale_color_manual(name = "ancestry \nmean block size", breaks = ancestryString, labels = ancestryString, values = ancestryColor) +
  labs(x = "generation", y = "average block size (bp)") +
  facet_wrap(~rec)
ggsave("averageLengthBlock.jpeg", width = 9*3, height = 9)








