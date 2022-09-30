## R version 4.1.1
## tidyverse version 1.3.1

## R script to plot the PCA results

library(tidyverse)

colorPopulation <- data.frame(col = c("#B8281D", "#7E8CBF", "#3C51A3", "#ABD39F", "#387E49", "#000001", "#C3C550", "#727374", "#BEBEC0", "#E598BC", "#D9366A"), pop = c("hpAfrica2","hspAfrica1SAfrica","hspAfrica1WAfrica","hspCNEAfrica","hspENEAfrica","hpAsia2","hspEAsia","hspEuropeC","hspEuropeN","hspEuropeSW","hspMiddleEast"))


pca <- read_table2("./coreCDS_LD_PCA.eigenvec", col_names = FALSE)
eigenval <- scan("./coreCDS_LD_PCA.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))


pca$ind <- sapply(strsplit(pca$ind, split = "_"), function(x) x[2])
pca$ind[1] <- "26695"
pop <- read.table("pylori_populations.tsv")

pca$ind[15] <- "26695"

pca$pop <- NA
for(i in 1:nrow(pca)) {
  pca$pop[i] <- pop$V2[pop$V1 == pca$ind[i]]
}

pop <- read.table("pop.txt", sep = "\t", fill = TRUE)
pop <- pop[!(pop$V4 %in% c("recipient","-")),]
pca$state <- "receiver"
pca$state[pca$ind %in% pop$V1] <- "donor"

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

colPop <- colorPopulation
pop <- unique(pca$pop)

# plot pca
fileName  <- paste0(folderPlots, "/PCA.jpeg")
p <- ggplot(pca, aes(PC1, PC2, col = pop, fill = pop, shape = state)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = c(20,3)) +
  scale_colour_manual(values = colorPopulation$col, breaks = colorPopulation$pop) +
  scale_fill_manual(values = colorPopulation$col, breaks = colorPopulation$pop) +
  xlab(paste0("PC1 (", round(pve$pve[1],1), "%)")) +
  ylab(paste0("PC2 (", round(pve$pve[2],1), "%)"))
p + theme_cowplot() + theme(text = element_text(size = 14),axis.text = element_text(size = 21),legend.justification = c(1, -0.1),legend.position = c(0.2, 0.3))
ggsave(fileName, width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)



pcal <- pca %>% pivot_longer(!c(ind,pop, state, PC1), names_to = "pc", values_to = "values")

# fileName  <- paste0(folderPlots, "/PCATot.jpeg")
# p <- ggplot(pcal, aes(PC1, values, col = pop, fill = pop, shape = state)) + 
#   geom_point(size = 3) +
#   scale_shape_manual(values = c(20,3)) +
#   scale_colour_manual(values = colorPopulation$col, breaks = colorPopulation$pop) +
#   scale_fill_manual(values = colorPopulation$col, breaks = colorPopulation$pop) +
#   facet_wrap(~ pc, scales = "free")
# p + theme_cowplot() + theme(text = element_text(size = 21),axis.text = element_text(size = 21),legend.justification = c(1, -0.1),legend.position = c(1, 0.4))
# ggsave(fileName, width = 8.5, height = 8.5, units = "in", dpi = 500, scale = 0.75)

fileName  <- paste0(folderPlots, "/PCA4.jpeg")
p <- ggplot(pca, aes(PC1, PC4, col = pop, fill = pop, shape = state)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = c(20,3)) +
  scale_colour_manual(values = colorPopulation$col, breaks = colorPopulation$pop) +
  scale_fill_manual(values = colorPopulation$col, breaks = colorPopulation$pop) +
  xlab(paste0("PC1 (", round(pve$pve[1],1), "%)")) +
  ylab(paste0("PC4 (", round(pve$pve[4],1), "%)"))
p + theme_cowplot() + theme(text = element_text(size = 14),axis.text = element_text(size = 21),legend.justification = c(1, -0.1),legend.position = c(0.2, 0.3))
ggsave(fileName, width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)

fileName  <- paste0(folderPlots, "/PCA3.jpeg")
p <- ggplot(pca, aes(PC1, PC3, col = pop, fill = pop, shape = state)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = c(20,3)) +
  scale_colour_manual(values = colorPopulation$col, breaks = colorPopulation$pop) +
  scale_fill_manual(values = colorPopulation$col, breaks = colorPopulation$pop) +
  xlab(paste0("PC1 (", round(pve$pve[1],1), "%)")) +
  ylab(paste0("PC3 (", round(pve$pve[3],1), "%)"))
p + theme_cowplot() + theme(text = element_text(size = 14),axis.text = element_text(size = 21),legend.justification = c(1, -0.1),legend.position = c(0.2, 0))
ggsave(fileName, width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)


