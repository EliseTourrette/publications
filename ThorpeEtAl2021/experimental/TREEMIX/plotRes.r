## script to plots the results of the treemix analysis
## note that I used hpAfrica2 as an outgroup

## R 4.1.1
## RColorBrewer 1.1-2
## R.utils 2.11.0
## OptM


library(RColorBrewer)
library(R.utils)
source("../../../plotting_funcs.R") 

prefix="core.noN"

for(edge in 0:15){
  pdf(paste0("../results_", edge,".pdf"))
  par(mfrow = c(2,1))
  for(rep in 1:10) {
    plot_tree(cex=0.8,paste0(prefix,".",edge,".",rep))
  }
  for(rep in 1:10) {
    plot_resid(stem=paste0(prefix,".",edge,".",rep),pop_order="../pop.list")
  }
}

graphics.off()


## inferring the optimal number of migration edges
## use the R package OptM?
library(OptM)

res <- optM(folder = ".", tsv = "outOptM")
plot_optM(res, pdf = "plotOptM")


## plot with the optimum number of edges
prefix="core.noN"
edge <- 9
rep <- 5

pdf(paste0(folderPlots, "/results_", edge,".pdf"), width = 16, height = 9)
plot_tree(cex=0.8,paste0(prefix,".",edge,".",rep))
graphics.off()


## look at the results of the threepop command that calculate the f3 statistics
## between pop A;B,C
## if f3 significantly less than 0 then pop A is admixed

data <- read.table("threepop.out", fill = TRUE)
colnames(data) <- c("pop", "f3", "se", "Z")
data <- data[-(1:2),]
data <- data[-which(data$pop == "Estimating" | data$pop == "total_nsnp"),1:4]

sub <- data[data$f3 < 0,]
sub <- sub[order(sub$f3),]



