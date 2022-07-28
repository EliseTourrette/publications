## r script to plot the fig showing the proportion of ancestry from each population
## separate plot for each recombination level
## add the proportion of unknown ancestry



setwd("/home/elise/Desktop/WORK/articleHP/simu10/simu10/analysesAncestry")

library(ggplot2)
library(cowplot)

selCoeff <- -c(0, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4)
bottleneck <- c(5000, 5500)
admixture <- seq(from = 8000, to = 15000, by = 500)
rec <- c('0bp', '500bp', '5000bp', '50000bp')


ancestryString <- c("bottleneck", "unknown", "non_bottleneck")
ancestryColor <- c("blue", "grey", "red")
ancestryLab <- c("bottlenecked \npopulation", "\nunknown \n", "non-bottlenecked \npopulation ")

load("proportion.RData")


irec <- "highRec"

p <- ggplot(proportionTot[!(proportionTot$population %in% c("pre_split")) & proportionTot$rec == irec,]) +
  geom_point(aes(x = G, y = proportion, color = population),alpha = 1,shape = 21,size = 1.5) +
  scale_color_manual(name = "ancestry", breaks = ancestryString, values = ancestryColor, labels = ancestryLab) +
  labs(x = "generation", y = "proportion") 
p + theme(text = element_text(size = 16),axis.text = element_text(size = 16),legend.justification = c(1, -0.1),legend.position = c(1, 0.35), panel.background = element_rect(fill = 'white', colour = 'white'), axis.line = element_line(size=0.5)) 
ggsave(paste0("proportion_", irec, ".jpeg"),width = 8.5,height = 8.5,units = "in",dpi = 500,scale = 0.75)






