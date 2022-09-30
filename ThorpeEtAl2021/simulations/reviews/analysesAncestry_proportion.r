## this script will analyse and plot the ancestry matrix 
## obtained by the python script ancestry.py (or ancestry2.py, 
## depending if want to look at only the polymorphic sites or at all sites)
## (because I am useless at handling matrix and plotting in python...)

## here calculate, the average proportion of the sites coming from each population
## do not discretize by when they arose

# setwd("/home/elise/Desktop/WORK/articleHP/simu10/simu10/analysesAncestry")
library(ggplot2)
library(cowplot)



ancestryString <- c("bottleneck", "unknown", "non_bottleneck")
ancestryColor <- c("blue", "grey", "red")

proportionTot <- NULL

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

proportion <- apply(data2+1, 1, tabulate, nbins = 6)
nbSite <- apply(proportion, 2, sum)
proportion <- data.frame(non_bottleneck = proportion[2,] + proportion[4,], bottleneck = proportion[3,] + proportion[5,], pre_split = proportion[6,], unknown = proportion[1,])
proportion <- 100 * proportion / nbSite

meanProp <- apply(proportion, 2, mean)
population <- attributes(meanProp)
attributes(meanProp) <- NULL
proportionTot <- rbind(proportionTot, data.frame(rec = rec, G = G, proportion = meanProp, population = population$names))

    
}


save(proportionTot, file = "proportion.RData")
#load("proportionData.RData")


ggplot(proportionTot[!(proportionTot$population %in% c("pre_split")),]) +
  geom_point(aes(x = G, y = proportion, color = population)) +
  scale_color_manual(name = "ancestry", breaks = ancestryString, values = ancestryColor) +
  labs(x = "generation", y = "proportion") +
  facet_wrap(~rec)
ggsave("proportion.jpeg", width = 9*3, height = 9)


ggplot(proportionTot[!(proportionTot$population %in% c("pre_split","unknown")),]) +
  geom_point(aes(x = G %% 500, y = proportion, color = population)) +
  scale_color_manual(name = "ancestry", breaks = ancestryString, values = ancestryColor) +
  labs(x = "generation", y = "proportion") +
  facet_wrap(~rec)
ggsave("proportion2.jpeg", width = 9*3, height = 9)


#admixture <- seq(from = 8000, to = 15000, by = 500)
#for(i in unique(proportionTot$rec)) {
#   fileName <- paste0(folderPlots, "/proportion", i, ".jpeg")
#    p <- ggplot(proportionTot[!(proportionTot$population %in% c("pre_split")) & (proportionTot$rec == i),]) +
#        geom_point(aes(x = G, y = proportion, color = population)) +
#        scale_color_manual(name = "ancestry", breaks = ancestryString, values = ancestryColor) +
#        labs(x = "generation", y = "ancestry proportion") +
#        geom_segment(data = data.frame(admixture = admixture), aes(x = admixture,y = -2.8,xend = admixture,yend = -0.5), arrow = arrow(length = unit(0.05, "in")))
#    p + theme_cowplot() + theme(text = element_text(size = 16),axis.text = element_text(size = 21),legend.justification = c(1, -0.1),legend.position = c(1, 0.4))
#    ggsave(fileName, width = 8.5, height = 8.5, units = "in", dpi = 500, scale = 0.75)
#}





