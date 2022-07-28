## all the figures for the articles and only the figures for the article

setwd("~/Desktop/articleRepeatedOUtOfAfricaPylori")

library(ggplot2)
library(cowplot)

selCoeff <- -c(0, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4)
bottleneck <- c(5000, 5500)
admixture <- seq(from = 8000, to = 15000, by = 500)
rec <- c('0bp', '500bp', '5000bp', '50000bp')

########### FIGURE 4 ###########

############# 4.B ###########

irec <- 3

load("~/Desktop/articleRepeatedOUtOfAfricaPylori/ANALYSIS/DATA/pairwiseDiff.RData")
data <- data[data$generation == 8000 & data$rec == rec[irec] & data$rep == 1, ]

selCoeff <- c(0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001, 0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001)
mutType <- paste0("m", 1:length(selCoeff))
for (i in 1:length(mutType)) {
  data$s[data$mutType == mutType[i]] <- selCoeff[i]
}

data1 <- data[data$mutType %in% mutType[1:7], -c(1:3, 5)]
data2 <- data[data$mutType %in% mutType[8:14], -c(1:3, 5)]

sub <- data.frame(pop = data1$pop, s = data1$s, data1[, 2:4951] + data2[, 2:4951])

data <- NULL
for (i in 2:7) {
  data <- rbind(data, data.frame(ratio = unlist(sub[i, 3:ncol(sub)]) / unlist(sub[1, 3:ncol(sub)]), neutral = unlist(sub[1, 3:ncol(sub)]), s = sub[i, 2], pop = sub[i, 1]))
  data <- rbind(data, data.frame(ratio = unlist(sub[i + 7, 3:ncol(sub)]) / unlist(sub[8, 3:ncol(sub)]), neutral = unlist(sub[8, 3:ncol(sub)]), s = sub[i + 7, 2], pop = sub[i + 7, 1]))
}

pathwd <- paste0("20201105_", irec, "_1")

data <- NULL
data <- rbind(data, data.frame(ratio = unlist(apply(sub[2:7, 3:ncol(sub)], 2, sum)) / unlist(sub[1, 3:ncol(sub)]), neutral = unlist(sub[1, 3:ncol(sub)]), pop = 1))
data <- rbind(data, data.frame(ratio = unlist(apply(sub[9:14, 3:ncol(sub)], 2, sum)) / unlist(sub[8, 3:ncol(sub)]), neutral = unlist(sub[8, 3:ncol(sub)]), pop = 2))

data2 <- aggregate(data[c("ratio", "neutral")], by = list(pop = data$pop), mean, na.rm = TRUE)

p <- ggplot(data) +
  geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)), alpha = 0.1, shape = 1, size = 0.75) +
  geom_point(data2,mapping = aes(x = neutral, y = ratio, fill = as.factor(pop)),alpha = 1,shape = 21,size = 1.5) +
  labs(x = "dS",y = "dN/dS",color = "Population",fill = "Population",linetype = " ") +
  scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation")) +
  scale_fill_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation"))
p + theme_cowplot() + theme(text = element_text(size = 21),axis.text = element_text(size = 21),legend.justification = c(1, -0.1),legend.position = c(1, 0.6))
ggsave(paste0("4B.jpg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)



############# 4.C ###########

i <- 3

pathwd <- paste0('20201105_', i, '_1')
load(file = paste0('~/Desktop/articleRepeatedOUtOfAfricaPylori/ANALYSIS/DATA/dNdS/dNdSancestorSumSelCoeff_', pathwd, '_g8000.RData'))

data2 <- aggregate(data[c("ratio", "neutral")], by = list(pop = data$pop), mean)

p <- ggplot(data) +
  geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)),alpha = 0.1,shape = 1,size = 0.75) +
  geom_point(data2,mapping = aes(x = neutral, y = ratio, fill = as.factor(pop)),alpha = 1,shape = 21,size = 1.5) +
  labs(x = "dS",y = "dN/dS",color = "Population",linetype = " ") +
  scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottleneck population", "bottleneck population")) +
  scale_fill_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottleneck population", "bottleneck population"))
p + theme_cowplot()  + theme(legend.position = "none",text = element_text(size = 21),axis.text = element_text(size = 21))
ggsave(paste0("4C.jpg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)




############# 4.D ###########

g <- 8300
irow <- 1
rep <- 1
cutOff <- 1000

irec <- 3

pathwd <- paste0("20201105_", irec, "_", rep)
load(paste0("/Users/elise/Desktop/articleRepeatedOUtOfAfricaPylori/ANALYSIS/DATA/ancestry_g15000/",pathwd,"/sampleData.RData"))

load(paste0("ANALYSIS/DATA/20211021-ancestry/dataAncestry",cutOff,"_",pathwd,".RData"))
anc <- data

anc[4:5, -(1:2)] <- anc[4:5, -(1:2)] / 10
## to get the ancestry in percent
## in generation 15000, corresponding to lines 5 to 7, the sample size is 1000

anc <- anc[, -(1:2)]

pos <- unlist(anc[irow + 2, ])
x1 <- unlist(anc[irow, ])
x2 <- unlist(anc[irow + 1, ])

data <- cbind(sub,ancestry1 = unlist(lapply(sub$pos, function(x) x1[pos == x])),ancestry2 = unlist(lapply(sub$pos, function(x) x2[pos == x])),pos2 = unlist(lapply(sub$pos, function(x) pos[pos == x])))

data$ratio1 <- data$ancestry1 / (data$ancestry1 + data$ancestry2)

## average over the delta fr bins
bins = c(-1.01,-0.99,-0.75,-0.5,-0.25,-0.01,0.01,0.25,0.5,0.75,0.99,1.01)
#bins = seq(from = -1.01, to = 1.01, by = 0.1)
data$bins <- NA
for (i in 2:length(bins)) {
  data$bins[data$fr > bins[i - 1] &data$fr <= bins[i]] <- mean(c(bins[i - 1], bins[i]))
}

p <-ggplot(data, aes(x = -bins,y = 1 - ratio1,col = as.factor(mtype))) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se,fun.args = list(mult = 1),geom = "errorbar",width = 0.05) +
  geom_point(mapping = aes(x = -fr,y = 1 - ratio1,col = factor(mtype)),alpha = 0.1) +expand_limits(y = c(0, 1.1)) +
  labs(x = "Mutation score", y = " Bottleneck population ancestry", color = "Mutation type") +
  geom_segment(aes(x = -0.01,y = 1.2,xend = -0.1,yend = 1.2),
  col = "black",
  arrow = arrow(length = unit(0.1, "in"))) +
  geom_text(x=-0.6, y=1.2, label="Excess \nbottleneck mutations", col = "black", size = 5) +
  geom_segment(aes(x = 0.01,y = 1.2,xend = 0.1,yend = 1.2),
  col = "black",
  arrow = arrow(length = unit(0.1, "in"))) +
  geom_text(x=0.58, y=1.2, label="Excess \nnon-bottleneck mutations", col = "black", size = 5)
p + theme_cowplot() + theme(text = element_text(size = 21),axis.text = element_text(size = 21))
ggsave(paste0("4D.jpg"),width = 10,height = 8.5,units = "in",dpi = 500,scale = 0.75)




############# FIGURE 4.A ###########

load("./ANALYSIS/DATA/fitnessOverall.RData")
colnames(data) <- c("rec", "rep", "pop", "generation", "mean", "var")

data$varnorm <- sqrt(data$var) / data$mean

#data <- data[-which(data$rec == "8e+05bp"),]

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

p <- ggplot(data[data$rec == unique(data$rec)[3],]) +
  stat_summary(mapping = aes(x = generation, y = mean, col = as.factor(pop)),fun.data = mean_se,geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  #facet_wrap(~ rec) +
  labs(x = "generation", y = "fitness", color = "Population") +
  scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation")) +
  geom_vline(xintercept = bottleneck[1], lty = "dotted") +
  geom_vline(xintercept = bottleneck[2], lty = "dotted")+
  geom_segment(data = data.frame(admixture = admixture), aes(x = admixture,y = 0.1 + seq(from = 0.15, to = 0,length.out = length(admixture)),xend = admixture,yend = 0.15 + seq(from = 0.15, to = 0,length.out = length(admixture))),
  arrow = arrow(length = unit(0.1, "in")))
p + theme_cowplot() + theme(text = element_text(size = 21),axis.text = element_text(size = 21),legend.justification = c(1, 0),legend.position = c(0.9, 0.65))
ggsave("4.A.jpg",width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)

############# FIGURE S6 ###########

############# S6.A ###########

irec <- 1

load("~/Desktop/articleRepeatedOUtOfAfricaPylori/ANALYSIS/DATA/pairwiseDiff.RData")
data <- data[data$generation == 8000 &data$rec == rec[irec] & data$rep == 1, ]

selCoeff <- c(0,-0.005,-0.002,-0.001,-0.0005,-0.0002,-0.0001,0,-0.005, -0.002,-0.001,-0.0005,-0.0002,-0.0001)
mutType <- paste0("m", 1:length(selCoeff))
for (i in 1:length(mutType)) {
  data$s[data$mutType == mutType[i]] <- selCoeff[i]
}

data1 <- data[data$mutType %in% mutType[1:7], -c(1:3, 5)]
data2 <- data[data$mutType %in% mutType[8:14], -c(1:3, 5)]

sub <-data.frame(pop = data1$pop, s = data1$s, data1[, 2:4951] + data2[, 2:4951])

data <- NULL
for (i in 2:7) {
  data <-rbind(data,data.frame(ratio = unlist(sub[i, 3:ncol(sub)]) / unlist(sub[1, 3:ncol(sub)]),neutral = unlist(sub[1, 3:ncol(sub)]),s = sub[i, 2],pop = sub[i, 1]))
  data <-rbind(data,data.frame(ratio = unlist(sub[i + 7, 3:ncol(sub)]) / unlist(sub[8, 3:ncol(sub)]),neutral = unlist(sub[8, 3:ncol(sub)]), s = sub[i + 7, 2], pop = sub[i + 7, 1]))
}

pathwd <- paste0("20201105_", irec, "_1")

data <- NULL
data <-rbind(data,data.frame(ratio = unlist(apply(sub[2:7, 3:ncol(sub)], 2, sum)) / unlist(sub[1, 3:ncol(sub)]), neutral = unlist(sub[1, 3:ncol(sub)]), pop = 1))
data <-rbind(data,data.frame( ratio = unlist(apply(sub[9:14, 3:ncol(sub)], 2, sum)) / unlist(sub[8, 3:ncol(sub)]), neutral = unlist(sub[8, 3:ncol(sub)]), pop = 2 ))

data2 <-aggregate(data[c("ratio", "neutral")], by = list(pop = data$pop), mean, na.rm = TRUE)

p <- ggplot(data) +
  geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)),alpha = 0.1, shape = 1, size = 0.75) +
  geom_point( data2,mapping = aes(x = neutral, y = ratio, fill = as.factor(pop)),alpha = 1, shape = 21,size = 1.5) +
  labs(x = "dS",y = "dN/dS",color = "Population",fill = "Population",linetype = " ") +
  scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c( "non-bottlenecked \npopulation \n", "bottlenecked \npopulation")) +
  scale_fill_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation"))
p + theme_cowplot() + theme(text = element_text(size = 21),axis.text = element_text(size = 21),legend.justification = c(1, -0.25),legend.position = c(0.9, 0.5))
ggsave(paste0("S6A.jpg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)




############# S6.B ###########

irec <- 4

load("~/Desktop/articleRepeatedOUtOfAfricaPylori/ANALYSIS/DATA/pairwiseDiff.RData")
data <-data[data$generation == 8000 &data$rec == rec[irec] & data$rep == 1, ]

selCoeff <- c(0,-0.005,-0.002,-0.001,-0.0005,-0.0002,-0.0001,0,-0.005, -0.002,-0.001,-0.0005,-0.0002,-0.0001)
mutType <- paste0("m", 1:length(selCoeff))
for (i in 1:length(mutType)) {
  data$s[data$mutType == mutType[i]] <- selCoeff[i]
}

data1 <- data[data$mutType %in% mutType[1:7], -c(1:3, 5)]
data2 <- data[data$mutType %in% mutType[8:14], -c(1:3, 5)]

sub <-data.frame(pop = data1$pop, s = data1$s, data1[, 2:4951] + data2[, 2:4951])

data <- NULL
for (i in 2:7) {
  data <-rbind(data,data.frame(ratio = unlist(sub[i, 3:ncol(sub)]) / unlist(sub[1, 3:ncol(sub)]), neutral = unlist(sub[1, 3:ncol(sub)]), s = sub[i, 2], pop = sub[i, 1]))
  data <-rbind(data,data.frame( ratio = unlist(sub[i + 7, 3:ncol(sub)]) / unlist(sub[8, 3:ncol(sub)]), neutral = unlist(sub[8, 3:ncol(sub)]), s = sub[i + 7, 2], pop = sub[i + 7, 1]))
}

pathwd <- paste0("20201105_", irec, "_1")

data <- NULL
data <-rbind(data,data.frame(ratio = unlist(apply(sub[2:7, 3:ncol(sub)], 2, sum)) / unlist(sub[1, 3:ncol(sub)]),neutral = unlist(sub[1, 3:ncol(sub)]),pop = 1))
data <-rbind(data,data.frame(ratio = unlist(apply(sub[9:14, 3:ncol(sub)], 2, sum)) / unlist(sub[8, 3:ncol(sub)]),neutral = unlist(sub[8, 3:ncol(sub)]),pop = 2))

data2 <-aggregate(data[c("ratio", "neutral")], by = list(pop = data$pop), mean, na.rm = TRUE)

p <- ggplot(data) +
  geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)),alpha = 0.1,shape = 1,size = 0.75) +
  geom_point(data2, mapping = aes(x = neutral, y = ratio, fill = as.factor(pop)),alpha = 1,shape = 21,size = 1.5) +
  labs(x = "dS",y = "dN/dS",color = "Population",linetype = " ") +
  scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation")) +
  scale_fill_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation"))
p + theme_cowplot() + theme(legend.position = "none",text = element_text(size = 21),axis.text = element_text(size = 21))
ggsave(paste0("S6B.jpg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)



############# S6.C ###########


i <- 1

pathwd <- paste0('20201105_', i, '_1')
load(file = paste0('~/Desktop/articleRepeatedOUtOfAfricaPylori/ANALYSIS/DATA/dNdS/dNdSancestorSumSelCoeff_', pathwd, '_g8000.RData'))

data2 <-
  aggregate(data[c("ratio", "neutral")], by = list(pop = data$pop), mean)

p <- ggplot(data) +
  geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)),alpha = 0.1,shape = 1,size = 0.75) +
  geom_point(data2,mapping = aes(x = neutral, y = ratio, fill = as.factor(pop)),alpha = 1,shape = 21,size = 1.5) +
  labs(x = "dS",y = "dN/dS",color = "Population",linetype = " ") +
  scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottleneck population", "bottleneck population")) +
  scale_fill_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottleneck population", "bottleneck population"))
p + theme_cowplot()  + theme(legend.position = "none",text = element_text(size = 21),axis.text = element_text(size = 21))
ggsave(paste0("S6C.jpg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)





############# S6.D ###########


i <- 4

pathwd <- paste0('20201105_', i, '_1')
load(file = paste0('~/Desktop/articleRepeatedOUtOfAfricaPylori/ANALYSIS/DATA/dNdS/dNdSancestorSumSelCoeff_', pathwd, '_g8000.RData'))

data2 <-aggregate(data[c("ratio", "neutral")], by = list(pop = data$pop), mean)

p <- ggplot(data) +
  geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)),alpha = 0.1,shape = 1,size = 0.75) +
  geom_point(data2,mapping = aes(x = neutral, y = ratio, fill = as.factor(pop)),alpha = 1,shape = 21,size = 1.5) +
  labs(x = "dS",y = "dN/dS",color = "Population",linetype = " ") +
  scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottleneck population", "bottleneck population")) +
  scale_fill_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottleneck population", "bottleneck population"))
p + theme_cowplot()  + theme(legend.position = "none",text = element_text(size = 21),axis.text = element_text(size = 21))
ggsave(paste0("S6D.jpg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)





############# S6.E ###########


g <- 8300
irow <- 1
rep <- 1
cutOff <- 1000

irec <- 1

pathwd <- paste0("20201105_", irec, "_", rep)
load(paste0("/Users/elise/Desktop/articleRepeatedOUtOfAfricaPylori/ANALYSIS/DATA/ancestry_g15000/",pathwd,"/sampleData.RData"))

load(paste0("ANALYSIS/DATA/20211021-ancestry/dataAncestry",cutOff,"_",pathwd,".RData"))
anc <- data

anc[4:5, -(1:2)] <- anc[4:5, -(1:2)] / 10
## to get the ancestry in percent
## in generation 15000, corresponding to lines 5 to 7, the sample size is 1000

anc <- anc[, -(1:2)]

pos <- unlist(anc[irow + 2, ])
x1 <- unlist(anc[irow, ])
x2 <- unlist(anc[irow + 1, ])

data <- cbind(sub,ancestry1 = unlist(lapply(sub$pos, function(x) x1[pos == x])),ancestry2 = unlist(lapply(sub$pos, function(x) x2[pos == x])),pos2 = unlist(lapply(sub$pos, function(x) pos[pos == x])))

data$ratio1 <- data$ancestry1 / (data$ancestry1 + data$ancestry2)

## average over the delta fr bins
bins = c(-1.01,-0.99,-0.75,-0.5,-0.25,-0.01,0.01,0.25,0.5, 0.75,0.99,1.01)
#bins = seq(from = -1.01, to = 1.01, by = 0.1)
data$bins <- NA
for (i in 2:length(bins)) {
  data$bins[data$fr > bins[i - 1] &data$fr <= bins[i]] <- mean(c(bins[i - 1], bins[i]))
}

p <- ggplot(data, aes(x = -bins, y = 1 - ratio1,col = as.factor(mtype))) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se,fun.args = list(mult = 1),geom = "errorbar",width = 0.05) +
  geom_point(mapping = aes(x = -fr,y = 1 - ratio1,col = factor(mtype)),alpha = 0.1) +
  expand_limits(y = c(0.9, 1.1)) +
  labs(x = "Mutation score", y = " Bottleneck population ancestry", color = "Mutation type")+
  geom_segment(aes(x = -0.01,y = 1.2,xend = -0.1,yend = 1.2),col = "black",arrow = arrow(length = unit(0.1, "in"))) +
  geom_text(x=-0.6, y=1.2, label="Excess \nbottleneck mutations", col = "black", size = 5) +
  geom_segment(aes(x = 0.01,y = 1.2,xend = 0.1,yend = 1.2),col = "black",arrow = arrow(length = unit(0.1, "in"))) +
  geom_text(x=0.58, y=1.2, label="Excess \nnon-bottleneck mutations", col = "black", size = 5)
p + theme_cowplot() + theme(text = element_text(size = 21),axis.text = element_text(size = 21))
ggsave(paste0("S6E.jpg"),width = 11,height = 8.5,units = "in",dpi = 500,scale = 0.75)







############# S6.F ###########


g <- 8300
irow <- 1
rep <- 1
cutOff <- 1000

irec <- 4

pathwd <- paste0("20201105_", irec, "_", rep)
load(paste0("/Users/elise/Desktop/articleRepeatedOUtOfAfricaPylori/ANALYSIS/DATA/ancestry_g15000/",pathwd,"/sampleData.RData"))

load(paste0("ANALYSIS/DATA/20211021-ancestry/dataAncestry",cutOff,"_",pathwd,".RData"))
anc <- data

anc[4:5, -(1:2)] <- anc[4:5, -(1:2)] / 10
## to get the ancestry in percent
## in generation 15000, corresponding to lines 5 to 7, the sample size is 1000

anc <- anc[, -(1:2)]

pos <- unlist(anc[irow + 2, ])
x1 <- unlist(anc[irow, ])
x2 <- unlist(anc[irow + 1, ])

data <- cbind(sub,ancestry1 = unlist(lapply(sub$pos, function(x) x1[pos == x])),ancestry2 = unlist(lapply(sub$pos, function(x) x2[pos == x])),pos2 = unlist(lapply(sub$pos, function(x) pos[pos == x])))

data$ratio1 <- data$ancestry1 / (data$ancestry1 + data$ancestry2)

## average over the delta fr bins
bins = c(-1.01,-0.99,-0.75,-0.5,-0.25,-0.01,0.01,0.25,0.5, 0.75,0.99,1.01)
#bins = seq(from = -1.01, to = 1.01, by = 0.1)
data$bins <- NA
for (i in 2:length(bins)) {
  data$bins[data$fr > bins[i - 1] & data$fr <= bins[i]] <- mean(c(bins[i - 1], bins[i]))
}

p <-ggplot(data, aes(x = -bins,y = 1 - ratio1,col = as.factor(mtype))) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_se,fun.args = list(mult = 1),geom = "errorbar",width = 0.05) +
  geom_point(mapping = aes(x = -fr,y = 1 - ratio1,col = factor(mtype)),alpha = 0.1) +
  expand_limits(y = c(0.9, 1.2)) +
  labs(x = "Mutation score", y = " Bottleneck population ancestry", color = "Mutation type")+
  geom_segment(aes(x = -0.01,y = 1.2,xend = -0.1,yend = 1.2),col = "black",arrow = arrow(length = unit(0.1, "in"))) +
  geom_text(x=-0.6, y=1.2, label="Excess \nbottleneck mutations", col = "black", size = 5) +
  geom_segment(aes(x = 0.01, y = 1.2,xend = 0.1, yend = 1.2),col = "black", arrow = arrow(length = unit(0.1, "in"))) +
  geom_text(x=0.58, y=1.1, label="Excess \nnon-bottleneck mutations", col = "black", size = 5)
p + theme_cowplot() + theme(text = element_text(size = 21),axis.text = element_text(size = 21))
ggsave(paste0("S6F.jpg"),width = 11,height = 8.5,units = "in", dpi = 500,scale = 0.75)




############# S6.G ###########

load("./ANALYSIS/DATA/fitnessOverall.RData")
colnames(data) <- c("rec", "rep", "pop", "generation", "mean", "var")

data$varnorm <- sqrt(data$var) / data$mean

#data <- data[-which(data$rec == "8e+05bp"),]

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

p <- ggplot(data[data$rec == unique(data$rec)[1],]) +
  stat_summary(mapping = aes(x = generation, y = mean, col = as.factor(pop)),fun.data = mean_se,geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  #facet_wrap(~ rec) +
  labs(x = "generation", y = "fitness", color = "Population") +
  scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation")) +
  geom_vline(xintercept = bottleneck[1], lty = "dotted") +
  geom_vline(xintercept = bottleneck[2], lty = "dotted")+
  geom_segment(data = data.frame(admixture = admixture), aes(x = admixture,y = 0.2 - seq(from = 0, to = 0.1,length.out = length(admixture)),xend = admixture,yend = 0.15 - seq(from = 0, to = 0.1,length.out = length(admixture))),arrow = arrow(length = unit(0.1, "in")))
p + theme_cowplot() + theme(text = element_text(size = 21),axis.text = element_text(size = 21),legend.position = "none")
ggsave("S6G.jpg",width = 8.5,height = 8.5, units = "in",dpi = 500,scale = 0.75)


############# S6.H ###########

load("./ANALYSIS/DATA/fitnessOverall.RData")
colnames(data) <- c("rec", "rep", "pop", "generation", "mean", "var")

data$varnorm <- sqrt(data$var) / data$mean

#data <- data[-which(data$rec == "8e+05bp"),]

data$rec <-as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

p <- ggplot(data[data$rec == unique(data$rec)[4],]) +
  stat_summary(mapping = aes(x = generation, y = mean, col = as.factor(pop)),fun.data = mean_se,geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  #facet_wrap(~ rec) +
  labs(x = "generation", y = "fitness", color = "Population") +
  scale_colour_manual(values = c("forestgreen", "Darkgoldenrod1"),labels = c("non-bottlenecked \npopulation \n","bottlenecked \npopulation")) +
  geom_vline(xintercept = bottleneck[1], lty = "dotted") +
  geom_vline(xintercept = bottleneck[2], lty = "dotted")+
  geom_segment(data = data.frame(admixture = admixture), aes(x = admixture,y = 0.58 + seq(from = 0.03, to = 0,length.out = length(admixture)),xend = admixture,yend = 0.6 + seq(from = 0.03, to = 0,length.out = length(admixture))),arrow = arrow(length = unit(0.1, "in")))
p + theme_cowplot() + theme(text = element_text(size = 21),axis.text = element_text(size = 21),legend.position = "none")
ggsave("S6H.jpg",width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)
