## figures for the H. pylori paper
## ancestry, dN/dS within population and to the ancestor 

library(ggplot2)
library(cowplot)

################### dN/dS within population ###################

rec <- c('0bp','500bp','5000bp','50000bp')

for(irec in c(1,3,4)) {
  
  load("~/Desktop/articleRepeatedOUtOfAfricaPylori/simu10_ResultsEffectBottleneckAdmixture/DATA/pairwiseDiff.RData")
  data <- data[data$generation == 8000 & data$rec == rec[irec] & data$rep == 1,]
  
  selCoeff <- c(0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001,0, -0.005, -0.002, -0.001, -0.0005, -0.0002, -0.0001)
  mutType <- paste0("m", 1:length(selCoeff))
  for(i in 1:length(mutType)) {
    data$s[data$mutType == mutType[i]] <- selCoeff[i]
  }
  
  data1 <- data[data$mutType %in% mutType[1:7],-c(1:3,5)]
  data2 <- data[data$mutType %in% mutType[8:14],-c(1:3,5)]
  
  sub <- data.frame(pop = data1$pop, s = data1$s, data1[,2:4951] + data2[,2:4951])
  
  data <- NULL
  for(i in 2:7) {
    data <- rbind(data, data.frame(ratio = unlist(sub[i,3:ncol(sub)])/unlist(sub[1,3:ncol(sub)]), neutral = unlist(sub[1,3:ncol(sub)]), s = sub[i,2], pop = sub[i,1]))
    data <- rbind(data, data.frame(ratio = unlist(sub[i+7,3:ncol(sub)])/unlist(sub[8,3:ncol(sub)]), neutral = unlist(sub[8,3:ncol(sub)]), s = sub[i+7,2], pop = sub[i+7,1]))
  }
  
  pathwd <- paste0("20201105_",irec,"_1")
  
  data <- NULL
  data <- rbind(data, data.frame(ratio = unlist(apply(sub[2:7,3:ncol(sub)],2,sum))/unlist(sub[1,3:ncol(sub)]), neutral = unlist(sub[1,3:ncol(sub)]), pop = 1))
  data <- rbind(data, data.frame(ratio = unlist(apply(sub[9:14,3:ncol(sub)],2,sum))/unlist(sub[8,3:ncol(sub)]), neutral = unlist(sub[8,3:ncol(sub)]), pop = 2))
  
  data2 <- aggregate(data[c("ratio","neutral")], by = list(pop = data$pop), mean, na.rm = TRUE)
  
  p <- ggplot(data) +
    geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)), alpha = 0.025, shape = 1) +
    geom_point(data2, mapping = aes(x = neutral, y = ratio, col = as.factor(pop)), alpha = 1, shape = 19) +  
    labs(x = "dS", y = "dN/dS", color = "Population", linetype = " ") +
    scale_colour_manual(values = c("#008081", "#D19900"), labels = c("non-bottleneck population", "bottleneck population")) 
  p + theme_cowplot() + theme(legend.position = "none")  
  ggsave(paste0("dNdSwithinPop",pathwd,".jpg"), width = 8.5, height = 8.5, units = "in",dpi = 500, scale = 0.75)
  
}

################### dN/dS to the ancestor ###################

for(i in c(1,3,4)) {
  pathwd <- paste0('20201105_',i,'_1')
  load(file = paste0('dNdS/dNdSancestorSumSelCoeff_',pathwd,'_g8000.RData'))
  
  data2 <- aggregate(data[c("ratio","neutral")], by = list(pop = data$pop), mean)
  
  p <- ggplot(data) +
    geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)), alpha = 0.025, shape = 1) +
    geom_point(data2, mapping = aes(x = neutral, y = ratio, col = as.factor(pop)), alpha = 1, shape = 19) +  
    labs(x = "dS", y = "dN/dS", color = "Population", linetype = " ") +
    scale_colour_manual(values = c("#008081", "#D19900"), labels = c("non-bottleneck population", "bottleneck population"))
  p + theme_cowplot()  + theme(legend.position = "none") 
  ggsave(paste0("dNdSancestor",pathwd,".jpg"), width = 8.5, height = 8.5, units = "in",dpi = 500, scale = 0.75)
  
}

################### ancestry ###################

g <- 8300
irow <- 1
rep <- 1
cutOff <- 1000
for(rec in c(1,3,4)) {
  
  pathwd <- paste0("20201105_",rec,"_",rep)
  load(paste0("/Users/elise/Desktop/articleRepeatedOUtOfAfricaPylori/simu10_ResultsEffectBottleneckAdmixture/ancestry_g15000/",pathwd,"/sampleData.RData"))
  
  load(paste0("20211021-ancestry/dataAncestry",cutOff,"_",pathwd,".RData"))
  anc <- data
  
  anc[4:5,-(1:2)] <- anc[4:5,-(1:2)]/10 
  ## to get the ancestry in percent 
  ## in generation 15000, corresponding to lines 5 to 7, the sample size is 1000
  
  anc <- anc[,-(1:2)]
  
  pos <- unlist(anc[irow+2,])
  x1 <- unlist(anc[irow,])
  x2 <- unlist(anc[irow+1,])
  
  data <- cbind(sub, 
                ancestry1 = unlist(lapply(sub$pos, function(x) x1[pos == x])), 
                ancestry2 = unlist(lapply(sub$pos, function(x) x2[pos == x])), 
                pos2 = unlist(lapply(sub$pos, function(x) pos[pos == x])))
  
  data$ratio1 <- data$ancestry1/(data$ancestry1 + data$ancestry2)
  
  ## average over the delta fr bins
  bins = c(-1.01,-0.99,-0.75,-0.5,-0.25,-0.01,0.01,0.25,0.5,0.75,0.99,1.01)
  #bins = seq(from = -1.01, to = 1.01, by = 0.1)
  data$bins <- NA
  for(i in 2:length(bins)) {
    data$bins[data$fr > bins[i-1] & data$fr <= bins[i]] <- mean(c(bins[i-1],bins[i]))
  }
  
  p <- ggplot(data, aes(x = -bins, y = 1 - ratio1, col = as.factor(mtype))) +
    stat_summary(fun = mean, geom = "point") +
    stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", width = 0.05) +
    geom_point(mapping = aes(x = -fr, y = 1 - ratio1, col = factor(mtype)), alpha = 0.1) +
    expand_limits(y=c(0,1)) +
    labs(x="mutation score", y = " Bottleneck population ancestry", color = "mutation type")
  p + theme_cowplot()
  ggsave(paste0("ancestry_g",g,"cutoff",cutOff,'_',pathwd,".jpg"), width = 8.5, height = 8.5, units = "in",dpi = 500, scale = 0.75)
}
















