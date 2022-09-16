## R version 4.1.1
## ggplot2 version 3.3.5

library(ggplot2)


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

ggplot(res, aes(x = zone, y = as.numeric(omegaNA), col = pop)) +
  geom_point() 

