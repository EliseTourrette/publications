## some additional figures dN/dS

######## DN/DS TO THE ANCESTOR

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

library(ggplot2)

for(i in c(1,3,4)) {
  pathwd <- paste0('20201105_',i,'_1')
  load(file = paste0('dNdSancestor_',pathwd,'_g8000.RData'))
  ggplot(data, aes(x = neutral, y = ratio, col = factor(s))) +
    geom_point()
  ggsave(paste0("dNdStoTheAncestor1",pathwd,".jpg"))
  
  ggplot(data, aes(x = neutral, y = ratio, col = pop)) +
    geom_point()
  ggsave(paste0("dNdStoTheAncestor2",pathwd,".jpg"))
  
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
  
  data2 <- aggregate(data[c("ratio","neutral")], by = list(pop = data$pop), mean)
  
  ggplot(data) +
    geom_point(mapping = aes(x = neutral, y = ratio, col = pop), alpha = 0.01) +
    geom_point(data2, mapping = aes(x = neutral, y = ratio, col = pop), alpha = 1) 
  ggsave(paste0("dNdStoTheAncestor3",pathwd,".jpg"))
  
}



##### DN/DS WITHIN POP 



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
  
  ggplot(data, aes(x = neutral, y = ratio, col = as.factor(s))) +
    geom_point()
  ggsave(paste0("dNdSwithinPop1",pathwd,".jpg"))
  
  ggplot(data, aes(x = neutral, y = ratio, col = as.factor(pop))) +
    geom_point()
  ggsave(paste0("dNdSwithinPop2",pathwd,".jpg"))
  
  
  data <- NULL
  data <- rbind(data, data.frame(ratio = unlist(apply(sub[2:7,3:ncol(sub)],2,sum))/unlist(sub[1,3:ncol(sub)]), neutral = unlist(sub[1,3:ncol(sub)]), pop = 1))
  data <- rbind(data, data.frame(ratio = unlist(apply(sub[9:14,3:ncol(sub)],2,sum))/unlist(sub[8,3:ncol(sub)]), neutral = unlist(sub[8,3:ncol(sub)]), pop = 2))
  
  data2 <- aggregate(data[c("ratio","neutral")], by = list(pop = data$pop), mean, na.rm = TRUE)
  
  ggplot(data) +
    geom_point(mapping = aes(x = neutral, y = ratio, col = as.factor(pop)), alpha = 0.01) +
    geom_point(data2, mapping = aes(x = neutral, y = ratio, col = as.factor(pop)), alpha = 1) 
  ggsave(paste0("dNdSwithinPop3",pathwd,".jpg"))
  
}
  





