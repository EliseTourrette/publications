## R version 4.1.1
## ggplot2 version 3.3.5

library(ggplot2)
library(cowplot)


folder <- "div_data"
files = unique(unlist(lapply(strsplit(system(paste("ls", folder), intern = TRUE), split = ".folded"), function(x) x[1])))
region <- c("Africa", "Asia", rep("Africa",2), "Asia", "Africa", rep("Europe / Middle East",4))

## use the scaledBeta model (4th model)
## still print the estimation for the neutral model (1st model)
res <- NULL
for(i in 1:length(files)) {
  fileData <- paste0(folder, "/", files[i], ".folded.dofe_basic.out")
  system(paste0("tail -n 8 ", fileData, " > tmp.txt"))
  data <- read.table(file = "tmp.txt", fill = TRUE)
  # AIC <- 2*(data$lnL) +2*data$X.param
  res <- rbind(res, data.frame(pop = files[i], model = data[6,1], omegaNA = data[6,2], zone = region[i]))
}

system("rm tmp.txt")

fileName <- paste0(folderPlots, "/omegaNA.jpeg")
p <- ggplot(res, aes(x = zone, y = as.numeric(omegaNA), col = pop)) +
    geom_point() +
    scale_colour_manual(values = colorPopulation$col, breaks = colorPopulation$pop) +
    xlab("region") +
    ylab("omega_NA")
p + theme_cowplot() + theme(text = element_text(size = 16),axis.text = element_text(size = 16),legend.justification = c(1, -0.1),legend.position = c(1, 0.4))
ggsave(fileName, width = 8.5, height = 8.5, units = "in", dpi = 500, scale = 0.75)



