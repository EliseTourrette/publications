## R version 4.1.1
## tidyverse version 1.3.1

## R script to plot the PCA results


library(tidyverse)
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

# plot pca
ggplot(pca, aes(PC1, PC2, col = pop, fill = pop, shape = state)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = c(20,3)) +
  xlab(paste0("PC1 (", round(pve$pve[1],1), "%)")) +
  ylab(paste0("PC2 (", round(pve$pve[2],1), "%)"))
ggsave("PCA.pdf", width = 16, height = 9)


pcal <- pca %>% pivot_longer(!c(ind,pop, state, PC1), names_to = "pc", values_to = "values")

ggplot(pcal, aes(PC1, values, col = pop, fill = pop, shape = state)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = c(20,3)) +
  facet_wrap(~ pc, scales = "free")
ggsave("PCATot.pdf", width = 16, height = 9)


ggplot(pca, aes(PC1, PC4, col = pop, fill = pop, shape = state)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = c(20,3)) +
  xlab(paste0("PC1 (", round(pve$pve[1],1), "%)")) +
  ylab(paste0("PC4 (", round(pve$pve[4],1), "%)"))
ggsave("PCA4.pdf", width = 16, height = 9)

ggplot(pca, aes(PC1, PC3, col = pop, fill = pop, shape = state)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = c(20,3)) +
  xlab(paste0("PC1 (", round(pve$pve[1],1), "%)")) +
  ylab(paste0("PC3 (", round(pve$pve[3],1), "%)"))
ggsave("PCA3.pdf", width = 16, height = 9)


