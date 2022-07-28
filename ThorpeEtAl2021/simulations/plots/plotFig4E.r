## script for the plot of FIG4E
## for the paper repeated out of africa ... H pylori
setwd("~/Desktop/articleRepeatedOutOfAfricaPylori/painting")

for (ig in c(8000, 8100, 8200, 8300)) {
  
  print(ig)
  
  nameFile <- paste0("genome2_", ig, ".txt") ## name of the file you will look at (number = generation)
  parameters <- read.table("parameters.csv", stringsAsFactor = FALSE, header = TRUE, sep = ",")
  
  nrow <- as.numeric(parameters$nrow[parameters$nameFile == nameFile])## number of rows of the file; get it from parameters.csv; can also get it with the "system" function 
  ncol <- as.numeric(system(paste("awk '{print NF}'", nameFile, "| sort -nu | tail -n 1"), intern = TRUE)) ## maximum number of column in the file
  
  ## load the population data
  pop <- read.table(nameFile, stringsAsFactor = FALSE, skip = nrow - 100, col.names = 1:ncol, fill = TRUE)
  pop <- as.matrix(pop[,-(1:2)])
  ## pop: matrix, rows = strains, columns = mutations contained by each strains
  
  ## load the mutation data
  mut <- read.table(nameFile,stringsAsFactor = FALSE, skip = 2,nrows = nrow - 100 - 3, fill = TRUE )
  mut <- mut[order(mut$V4),]
  ## mut: data.frame, rows = mutations, column = different information about the given mutations
  
  ##if want to look at a portion of the genome (or the whole genome), independently on whether a site has a mutation or not
  nbSite <- 1600000 ## if want to look at the whole genome, nbSite = 1600000
  ancestry <- matrix(0, nrow = nrow(pop), ncol = nbSite) ## if want to look at a portion of the chromosomes of nbSite sites
  posBegin <- 1 ## first position, in bp of the interval looked at
  
  ##fill the ancestry matrix
  for (i in 1:nbSite) {
    ## for-loop on each site
    j <- i + posBegin - 1 ## get the physical position of each site = index of the site in the interval plus the first position of the interval (minus 1, else would add 1 to every position)
    id <- mut$V1[mut$V4 == j] ## get the (within-file) id of the mutation(s) at this position (there can be more than one or none at all)
    if (length(id) > 0) {
      ## do the following part only if there is at least one mutation at this position
      #print("mut!")
      for (k in id) {
        ## as there can be more than mutation at this position, loop on the possible mutation
        ind <- which(pop == k) %% 100 ## get the index of the strains with the mutation (use the modulo 100 as there are 100 strains in the sample; would need another way if the number of strains in the sample was different)
        ind[ind == 0] <- 100 ## for the element that are a multiple of 100, the modulo will result in 0; thus need to change the strain index back to 100 (could also do ind <- ind + 1, as we need an index that begin at 1)
        # fitness[ind,i] <- 1 + mut$V5[mut$V1 == k]
        ##fill in the ancestry matrix depending on the mutation
        if (mut$V7[mut$V1 == k] == "p1" & mut$V8[mut$V1 == k] < 5000)
          ancestry[ind, i] <- 1
        if (mut$V7[mut$V1 == k] == "p1" & mut$V8[mut$V1 == k] >= 5000 & mut$V8[mut$V1 == k] < 8000)
          ancestry[ind, i] <- 2
        if (mut$V7[mut$V1 == k] == "p1" & mut$V8[mut$V1 == k] >= 5000 & mut$V8[mut$V1 == k] >= 8000)
          ancestry[ind, i] <- 3
        if (mut$V7[mut$V1 == k] == "p2" & mut$V8[mut$V1 == k] >= 5000 & mut$V8[mut$V1 == k] < 8000)
          ancestry[ind, i] <- 4
        if (mut$V7[mut$V1 == k] == "p2" & mut$V8[mut$V1 == k] >= 5000 & mut$V8[mut$V1 == k] >= 8000)
          ancestry[ind, i] <- 5
      }
    }
  }
  positionSite <- posBegin:(posBegin + nbSite - 1) ## vector of position, in bp, of continuous sites
  
  ## initialize the selcoeff matrix with 0, which is the default for a site in a strain without mutation
  selcoeff <- matrix(0, nrow = nrow(pop), ncol = length(unique(mut$V4)))
  ## put the selective coefficient in selcoeff matrix, rows = strains and columns = site
  for (i in 1:length(unique(mut$V4))) {
    ## for site which has a mutation; for-loop on a index
    j <- unique(mut$V4)[i] ## get the physical position
    id <- mut$V1[mut$V4 == j] ## get the (within-file) id of the mutation(s) at this position (there can be more than one)
    for (k in id) {
      ## as there can be more than mutation at this position, loop on the possible mutation
      ind <- which(pop == k) %% 100 ## get the index of the strains with the mutation (use the modulo 100 as there are 100 strains in the sample; would need another way if the number of strains in the sample was different)
      ind[ind == 0] <- 100 ## for the element that are a multiple of 100, the modulo will result in 0; thus need to change the strain index back to 100 (could also do ind <- ind + 1, as we need an index that begin at 1)
      ## fill in the selcoeff matrix depending on the mutation type(mut$3)
      selcoeff[ind, i] <- mut$V5[mut$V1 == k]
    }
  }
  
  ## matrix of fitness; fitness=1+selective coefficient
  fitness <- selcoeff + 1
  
  ## if want to look at the fitness for the different classes of mutations
  # components <- apply(fitness, 1, table) ## first, for each strain, calculate the number of fitness of each class of mutations
  # effect <- as.numeric(attributes(components)$dimnames[[1]]) ## the get the fitness effect (1+s) of each class of mutations present
  # fitnessPerClass <- effect^components ## then calculate the fitness per strain per clas of mutations (which is simply effect^(nbr mutations of this class) as the fitness components are multiplied between the different sites)
  
  ##fitness_i = mult{from j = 1 to N} (1 + s_ij), with i = strain, j = site (N = number of site) and s_ij the selection coeff of strain i at site j
  fitnessTotal <- apply(fitness, 1, prod)## here, calculate the total fitness of each strain
  
  ## for the test, sample 10 strains among the 100 strains
  Nstrains <- c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91)
  
  ancestry10 <- ancestry[Nstrains,]
  fitness10 <- fitnessTotal[Nstrains]
  
  widthBlock <- apply(ancestry10, 2, sum) ## if the value > 0, then there is at least one mutation at this site for this subset of strain
  widthBlock[widthBlock > 0] <- 1 ## only use 0 and 1
  
  library(reshape2)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
  
  anc10 <- melt(ancestry10) ## need to reshape the matrix to get a data.frame with col1 = strain, col2 = site, col3 = ancestry (note that you could directly put it in this format instead of the ancestry matrix
  ## add another variable to anc10 saying if there is a mutation at the site or not
  anc10$Var3 <- rep(widthBlock, each = length(Nstrains)) ## use each as we first look at every strains for a position; would use times if we were looking first at every site for one strain
  ## replace the site number by their physical positions
  anc10$Var2 <- rep(positionSite, each = length(Nstrains)) ## use each as we first look at every strains for a position; would use times if we were looking first at every site for one strain
  
  ## want to look at the global proportion mutations non bottleneck/bottleneck population for the mutations post split / pre admixture
  origin10 <- NULL
  for (i in 1:length(Nstrains)) {
    origin <- anc10$value[anc10$Var1 == i]
    origin10 <- c(origin10, length(origin[origin == 4]) / (length(origin[origin == 2])+length(origin[origin == 4])))
  }
  print(max(origin10))
  
  anc10 <- anc10[anc10$Var2 <= 20000,] ## only look at the first 20000 sites
  
  lab <- c( "no mutations",
            "pop1 pre split \n",
            "pop1, post split \npre admixture \n",
            "pop1, post split \npost admixture \n",
            "pop2, post split \npre admixture \n",
            "pop2, post split \npost admixture \n") ## possible labels for the ancestry
  colLab <- c("white",  "green4",  "red", "orange", "blue", "purple") ## color to be used for the different labels; same color scheme
  
  plot1 <- ggplot(data = anc10, aes(x = Var2, y = as.factor(Var1))) +
    scale_fill_manual(labels = lab[sort(unique(anc10$value + 1))], values = colLab[sort(unique(anc10$value + 1))]) +
    geom_tile(aes(fill = as.factor(value), height = 0.7, width = Var3)) +  ## add the plus if want to add one of the geom_point ## with width = Var3, you only put a bar when there is a mutation
    # geom_raster(aes(fill = as.factor(value))) +
    #geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), pch = as.factor(round(value,4)))) ## for each element of the ancestry matrix, will add a point whose shape will depend on the value of the selection coefficient
    # geom_point(data = s10, aes(x = Var2, y = as.factor(Var1), size = as.factor(abs(round(value,4)))))  ## for each element of the ancestry matrix, will add a point whose size will depend on the value of the selection coefficient
    theme_set(theme_bw()) + theme(
      panel.grid = element_line(colour = NA),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      #axis.text.x = element_text(size = rel(0.5), angle = 45),
      axis.line.x = element_line(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank()) + ##get rid of the grey background and grid lines
    guides(fill = FALSE)
  
  plot2 <-
    ggplot(data = data.frame(strain = 1:length(fitness10), fitness = fitness10), aes(x = fitness, y = as.factor(strain))) +
    geom_point(size = 0.5) +
    theme_set(theme_bw()) + theme(
      panel.grid = element_line(colour = NA),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      #axis.text.x = element_text(size = rel(0.5), angle = 45),
      panel.border = element_blank(),
      axis.line = element_line()) + ##get rid of the grey background and grid lines
    xlim(0.31, 0.37)
  
  plot3 <- ggplot(data = data.frame(strain = 1:length(origin10), ratio = origin10), aes(x = ratio, y = as.factor(strain))) +
    geom_point(size = 0.5) +
    theme_set(theme_bw()) + theme(
      panel.grid = element_line(colour = NA),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      #axis.text.x = element_text(size = rel(0.5), angle = 45),
      panel.border = element_blank(),
      axis.line = element_line()) + ##get rid of the grey background and grid lines
    xlim(0,1)
    
  plotTot <- plot_grid( plot2, NULL,  plot3, NULL, plot1,  nrow = 1, rel_widths = c(1, -0.01, 1, -0.01, 4),  align = "hv")
  
  # now add the title
  title <- ggdraw() + draw_label(paste0("g = ", ig), size = 12)
  p <- plot_grid(title, plotTot,  ncol = 1, rel_heights = c(1, 9))
  
  # save_plot(paste0("generation8000.pdf"), p, base_height = 1.8)
  
  save_plot(paste0("generation", ig, ".pdf"), p, base_height = 2.5,  dpi = 2000)
  
}
