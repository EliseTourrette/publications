## plot from the data save using analysisBottleneck3.r

library(ggplot2)
library(cowplot)

#setwd("C:/Users/Elise/Desktop/POSTDOC_IPS/Simulations/Runs/bottleneck/allReps")

selCoeff <- -c(0,5e-3,2e-3,1e-3,5e-4,2e-4,1e-4)
bottleneck <- c(5000,5500)





################# fitness per mutation class #################

load("./fitnessPerClass.RData")

colnames(data) <- c("rec", "rep", "mutType", "selCoeff", "pop", "generation", "fitness", "simu")

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

## format data: # generation
##              fitness, selection coefficient (for the deleterious mutations)
##              population, recombination level, replicate

ggplot(data, aes(x = generation, y = fitness, col = as.factor(selCoeff), fill = as.factor(selCoeff), linetype = as.factor(pop))) +
  stat_summary(fun.data = mean_se, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  facet_wrap(~ rec) +
  labs(x = "generation", y = "fitness", fill = "selection coefficient", color = "selection coefficient", linetype = "population") +
  scale_linetype_discrete(labels = c("ancestral population", "bottleneck population", "admixed population")) + 
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0.35, yend = 0.35), colour = "black") 
ggsave("averageFitnessClasses.jpg", width = 16, height = 9, units = "in", scale = 0.5)


################# overall fitness (mean and variance) #################

load("indpdtVar.RData")
data0 <- data
colnames(data0) <- c("rec", "rep", "pop", "generation", "var")


data0$rec <- as.numeric(unlist(strsplit(as.character(data0$rec), split = "bp")))
data0$rec <- paste("import size per generation: ", data0$rec, "bp")
data0$pop <- paste0("p", data0$pop)

load("./fitnessOverall.RData")
colnames(data) <- c("rec", "rep", "pop", "generation", "mean", "var")

data$varnorm <- sqrt(data$var) / data$mean

#data <- data[-which(data$rec == "8e+05bp"),]

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

p <- ggplot(data[data$rec == unique(data$rec)[4], ]) +
  stat_summary(mapping = aes(x = generation, y = mean, col = as.factor(pop)), fun.data = mean_se, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  #facet_wrap(~ rec) +
  labs(x = "generation", y = "fitness", color = "population") +
  scale_colour_manual(values = c("#008081", "#D19900", "red"), labels = c("ancestral population", "bottleneck population", "admixed population")) +
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0, yend = 0), colour = "black") 
p + theme_cowplot()+ theme(legend.position = "none") 
ggsave("averageFitness4.jpg", width = 8.5, height = 8.5, units = "in",dpi = 500, scale = 0.75)

data$pop <- paste0("p", data$pop)

ggplot(data) +
  stat_summary(mapping = aes(x = generation, y = var, col = as.factor(pop), lty = as.factor(1)), fun.data = mean_se, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  facet_wrap(~ rec) +
  labs(x = "generation", y = "fitness variance", color = "population", linetype = " ") +
  scale_colour_manual(values = c("#008081", "#D19900", "red"), labels = c("ancestral population", "bottleneck population", "admixed population")) +
  stat_summary(data = data0, mapping = aes(x = generation, y = var, col = as.factor(pop), lty = as.factor(2)), fun.data = mean_se, geom = "line") +
  scale_linetype_discrete(labels = c("linked sites","independent sites")) +
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0, yend = 0), colour = "black")  
ggsave("varFitness.jpg", width = 16, height = 9, units = "in", scale = 0.5)



dataSub <- data[which(data$generation %in% unique(data0$generation)),]
index <- with(dataSub, order(rep,rec, pop, generation))
dataSub <- dataSub[index,]

index <- with(data0, order(rep,rec, pop, generation))
data0 <- data0[index,]
data0$varNorm <- sqrt(data0$var)/dataSub$mean

p <- ggplot(data) +
  stat_summary(mapping = aes(x = generation, y = varnorm, col = as.factor(pop), lty = as.factor(1)), fun.data = mean_se, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  facet_wrap(~ rec) +
  labs(x = "generation", y = "fitness sd / mean", color = "population", linetype = " ") +
  scale_colour_manual(values = c("#008081", "#D19900"), labels = c("ancestral population", "bottleneck population")) +
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0, yend = 0), colour = "black") +
  scale_linetype_discrete(labels = c("simulated genomes","randomized genomes")) +
  stat_summary(data = data0, mapping = aes(x = generation, y = varNorm, col = as.factor(pop), lty = as.factor(2)), fun.data = mean_se, geom = "line") 
p + theme_cowplot()
ggsave("sdFitnessNorm.jpg", width = 16, height = 9, units = "in", scale = 0.75)



##########################################
## dN/dS to the ancestor

load("diffAncestor.RData")

data$selCoeff <- as.numeric(data$selCoeff)

#data <- data[-which(data$rec == "8e+05bp"),]
neutral <- data[which(data$selCoeff == 0),]
data <- data[ - which(data$selCoeff == 0),]

index <- with(data, order(selCoeff, rep,rec, pop, generation))
data <- data[index,]

index <- with(neutral, order(selCoeff, rep,rec, pop, generation))
neutral <- neutral[index,]

data$selCoeff <- - data$selCoeff
data$neutral <- neutral$diffAncestor
data$ratio <- data$diffAncestor/data$neutral * 6

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

p <- ggplot(data, aes(x = generation, y = ratio, col = as.factor(selCoeff), linetype = as.factor(pop))) +
  stat_summary(fun.data = mean_se, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  labs(x = "generation", y = "dN/dS to the ancestor", color = "selection coefficient", linetype = "population") +
  scale_linetype_discrete(labels = c("ancestral population", "bottleneck population", "admixed population")) +
  facet_wrap(~ rec) +
  geom_hline(yintercept = 1, lty = 2, colour = "grey") +
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0, yend = 0), colour = "black") 
p + theme_cowplot()
ggsave(paste("dNdSAncestor.jpg", sep = ""), width = 16, height = 9, units = "in", scale = 0.75)


datag <- data[data$generation == 15000 & data$rec == "import size per generation:  5000 bp",]
datag <- aggregate(datag[c("ratio","neutral")], by = list(selCoeff = datag$selCoeff, pop = datag$pop), mean)
p <- ggplot(datag, aes(x = neutral, y = ratio, col = as.factor(selCoeff), pch = as.factor(pop))) +
  stat_summary(fun.data = mean_se, geom = "point") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  labs(x = "dS", y = "dN/dS to the ancestor", color = "selection coefficient", linetype = "population") +
  scale_shape_discrete(labels = c("ancestral population", "bottleneck population", "admixed population"))  
p + theme_cowplot()
ggsave(paste("dNdSAncestor_dS.jpg", sep = ""), width = 16, height = 9, units = "in", scale = 0.75)


# datag <- data[data$generation == 10000,]
# 
# # datag$neutral <- datag$neutral/1.6e6 ## if want the number of differences compared to the genome size
# tdata <- aggregate(datag[c("diffAncestor", "neutral", "ratio")], by = list(selCoeff = datag$selCoeff, pop = datag$pop, generation = datag$generation, rec = datag$rec), mean, na.rm = TRUE)
# tdatag <- aggregate(datag[c("diffAncestor", "neutral", "ratio")], by = list(selCoeff = datag$selCoeff, pop = datag$pop, generation = datag$generation, rec = datag$rec), sd, na.rm = TRUE)
# 
# ggplot(tdata, aes(x = neutral/1.6e6, y = ratio, col = as.factor(selCoeff), pch = as.factor(pop))) +
#   stat_summary(fun = mean, geom = "point") +
#   #geom_errorbar(aes(xmin = neutral - 1.96*tdatag$neutral/2, xmax = neutral + 1.96*tdatag$neutral/2, y = ratio), width = 0.01) +
#   #geom_errorbar(aes(x = neutral, ymin = ratio - 1.96*tdatag$ratio/2, ymax = ratio + 1.96*tdatag$neutral/2), width = 0.01) +
#   labs(x = "dS", y = "dN/dS to the ancestor", title = "g = 3000; all classes of mutations", color = "selection coefficient", pch = "population") +
#   scale_shape_discrete(labels = c("ancestral population", "bottleneck population", "admixed population")) +
#   geom_hline(yintercept = 1, lty = 2, colour = "grey") +
#   facet_wrap(~ rec)
# ggsave(paste("dNdSAncestor2_allClasses.pdf", sep = ""), width = 16, height = 9, units = "in", scale = 0.5)
# 


##########################################
## dN/dS within pop

load("pairwiseDiff.RData")
colnames(data) <- c("rec", "rep", "mutType", "pop", "g", 1:4950)

#data <- data[-which(data$rec == "8e+05bp"),]

selCoeff <- c(0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001,0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001)
mutType <- paste0("m", 1:length(selCoeff))
for(i in 1:length(mutType)) {
  data$s[data$mutType == mutType[i]] <- selCoeff[i]
}


data0 <- data

data <- aggregate(data[,6:4955], by = list(rec = data$rec, rep = data$rep, pop = data$pop, g = data$g, s = data$s), sum)

data1 <- apply(data[,-c(1:5)], 1, mean, na.rm = TRUE)
data <- data[,c(1:5)]
data$m <- data1

neutral <- data[which(data$s == 0),]
data <- data[-which(data$s == 0),]

index <- with(data, order(s, rep,rec, pop, g))
data <- data[index,]

index <- with(neutral, order(s, rep,rec, pop, g))
neutral <- neutral[index,]

data$neutral <- neutral$m
data$ratio <- data$m/data$neutral * 6

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

p <- ggplot(data, aes(x = g, y = ratio, col = as.factor(s), linetype = as.factor(pop))) +
  stat_summary(fun.data = mean_se, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  labs(x = "generation", y = "dN/dS within pop", color = "selection coefficient", linetype = "population") +
  scale_linetype_discrete(labels = c("ancestral population", "bottleneck population", "admixed population")) +
  geom_hline(yintercept = 1, lty = 2, colour = "grey") +
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0, yend = 0), colour = "black") +
  facet_wrap(~ rec)
p + theme_cowplot()
ggsave(paste("dNdSwithinPop.jpg", sep = ""), width = 16, height = 9, units = "in", scale = 0.75)

data <- data0

data <- data[data$pop == 2,]

data$pop[data$mutType %in% paste0("m",1:7)] <- 1
data$pop[data$mutType %in% paste0("m",8:14)] <- 2

data1 <- apply(data[,-c(1:5,4956)], 1, mean, na.rm = TRUE)
data <- data[,c(1:5,4956)]
data$moy <- data1

neutral <- data[which(data$s == 0),]
data <- data[-which(data$s == 0),]

n1 <- neutral[neutral$pop == 1,]
n2 <- neutral[neutral$pop == 2,]
d1 <- data[data$pop == 1,]
d2 <- data[data$pop == 2,]

index <- with(d1, order(s, rep,rec, pop, g))
d1 <- d1[index,]
index <- with(d2, order(s, rep,rec, pop, g))
d2 <- d2[index,]

index <- with(n1, order(s, rep,rec, pop, g))
n1 <- n1[index,]
index <- with(n2, order(s, rep,rec, pop, g))
n2 <- n2[index,]

d1$neutral <- n1$moy
d2$neutral <- n2$moy
data <- rbind(d1,d2)
data$ratio <- data$moy/data$neutral * 6

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

ggplot(data, aes(x = g, y = ratio, col = as.factor(s), linetype = as.factor(pop))) +
  stat_summary(fun.data = mean_se, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  labs(x = "generation", y = "dN/dS within pop", color = "selection coefficient", linetype = "population") +
  geom_hline(yintercept = 1, lty = 2, colour = "grey") +
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0, yend = 0), colour = "black") +
  facet_wrap(~ rec, scales = "free")
ggsave(paste("dNdSwithinPop_pop2.jpg", sep = ""), width = 16, height = 9, units = "in", scale = 0.5)


datag <- data[data$generation == 15000 & data$rec == "import size per generation:  5000 bp",]
datag <- aggregate(datag[c("ratio","neutral")], by = list(selCoeff = datag$selCoeff, pop = datag$pop), mean)
p <- ggplot(datag, aes(x = neutral, y = ratio, col = as.factor(selCoeff), pch = as.factor(pop))) +
  stat_summary(fun.data = mean_se, geom = "point") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  labs(x = "dS", y = "dN/dS to the ancestor", color = "selection coefficient", linetype = "population") +
  scale_shape_discrete(labels = c("ancestral population", "bottleneck population", "admixed population"))  
p + theme_cowplot()
ggsave(paste("dNdSAncestor_dS.jpg", sep = ""), width = 16, height = 9, units = "in", scale = 0.75)


# datag <- data[data$g == max(data$g),]
# 
# tdata <- aggregate(datag[c("moy", "neutral", "ratio")], by = list(s = datag$s, pop = datag$pop, g = datag$g, rec = datag$rec), mean, na.rm = TRUE)
# tdatasd <- aggregate(datag[c("moy", "neutral", "ratio")], by = list(s = datag$s, pop = datag$pop, g = datag$g, rec = datag$rec), sd, na.rm = TRUE)
# 
# ggplot(tdata, aes(x = neutral/1.6e6, y = ratio, col = as.factor(s), pch = as.factor(pop))) +
#   stat_summary(fun = mean, geom = "point") +
#   labs(x = "dS", y = "dN/dS within pop", title = "g = 3000, all classes of mutations", color = "selection coefficient", pch = "population") +
#   scale_shape_discrete(labels = c("ancestral population", "bottleneck population", "admixed population")) +
#   geom_hline(yintercept = 1, lty = 2, colour = "grey") +
#   facet_wrap(~ rec)
# ggsave(paste("dNdSwithinPop2_allClasses.pdf", sep = ""), width = 16, height = 9, units = "in", scale = 0.5)
# 

##########################################
## pairwise differences

load("pairwiseDiff.RData")
colnames(data) <- c("rec", "rep", "mutType", "pop", "g", 1:4950)

selCoeff <- c(0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001,0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001)
mutType <- paste0("m", 1:length(selCoeff))
for(i in 1:length(mutType)) {
  data$s[data$mutType == mutType[i]] <- selCoeff[i]
}

data0 <- data

data <- aggregate(data[6:4955], by = list(rec = data$rec, rep = data$rep, pop = data$pop, g = data$g, s = data$s), sum)

data1 <- apply(data[,-c(1:5)], 1, mean, na.rm = TRUE)
data <- data[,c(1:5)]
data$m <- data1

#data <- data[-which(data$rec == "8e+05bp"),]
data$m[data$s == 0] <- data$m[data$s == 0]/6

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")


ggplot(data, aes(x = g, y = m, col = as.factor(s), linetype = as.factor(pop))) +
  stat_summary(fun = mean, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  facet_wrap(~ rec) +
  labs(x = "generation", y = "average pairwise difference", color = "selection coefficient", linetype = "population") +
  scale_linetype_discrete(labels = c("ancestral population", "bottleneck population", "admixed population")) +
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0, yend = 0), colour = "black")
ggsave("pairwiseDiff.jpg", width = 16, height = 9, units = "in", scale = 0.5)

data <- data0

data1 <- apply(data[,-c(1:5,4956)], 1, mean, na.rm = TRUE)
data <- data[,c(1:5,4956)]
data$m <- data1

#data <- data[-which(data$rec == "8e+05bp"),]
data$m[data$s == 0] <- data$m[data$s == 0]/6

data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

data <- data[data$pop == 2,]

data$pop[data$mutType %in% paste0("m",1:7)] <- 1
data$pop[data$mutType %in% paste0("m",8:14)] <- 2

ggplot(data, aes(x = g, y = m, col = as.factor(s), linetype = as.factor(pop))) +
  stat_summary(fun.data = mean_se, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  facet_wrap(~ rec) +
  labs(x = "generation", y = "average pairwise difference", color = "selection coefficient", linetype = "population of origin") +
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0, yend = 0), colour = "black")
ggsave("pairwiseDiff_pop2.jpg", width = 16, height = 9, units = "in", scale = 0.5)

##########################################
## differences with the ancestor: number of mutations PER strains

load("diffAncestor.RData")

data$selCoeff <- as.numeric(data$selCoeff)

#data <- data[-which(data$rec == "8e+05bp"),]
data$diffAncestor[data$selCoeff == 0] <- data$diffAncestor[data$selCoeff == 0]/6


data$rec <- as.numeric(unlist(strsplit(as.character(data$rec), split = "bp")))
data$rec <- paste("import size per generation: ", data$rec, "bp")

ggplot(data, aes(x = generation, y = diffAncestor, col = as.factor(selCoeff), linetype = as.factor(pop), fill = as.factor(selCoeff))) +
  stat_summary(fun.data = mean_se, geom = "line") +
  #stat_summary(fun.data = mean_se, fun.args = list(mult = 1.96), geom = "ribbon", alpha = 0.2) +
  facet_wrap(~ rec) +
  labs(x = "generation", y = "nbr of mutation per strains", fill = "selection coefficient", color = "selection coefficient", linetype = "population") +
  scale_linetype_discrete(labels = c("ancestral population", "bottleneck population", "admixed population")) +
  geom_segment(aes(x = bottleneck[1], xend = bottleneck[2], y = 0, yend = 0), colour = "black") 
ggsave("diffAncestor.jpg", width = 16, height = 9, units = "in", scale = 0.5)



##########################################
## frequencies of the mutations introduced by migration

load("timeFr.RData")
dm <- aggregate(data["fr"], by = list(rec = data$rec, g = data$g, s = data$s, gm = data$gm), mean) ## last generation of introduction into pop2 of the mutation by migration/bottleneck

rec = "50000bp"

ggplot(dm[dm$rec == rec,]) + 
  geom_point(mapping = aes(x = g, y = fr, col = as.factor(s)), shape= 20)+
  facet_wrap(~gm)
ggsave(paste0("plot_", rec, ".jpg"), width = 16, height = 9, units = "in", scale = 0.5) 





