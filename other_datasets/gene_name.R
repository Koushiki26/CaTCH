#gencode
setwd("/home/dell/Documents/koushiki/project_new/gencode/")

data1 <- read.csv("data6.csv", header=T)
head(data1)
data2 <- read.csv("gene_name.txt", header=F, sep = "\t")
head(data2)
colnames(data2) <- c("Gene_Id", "Gene_name")
library(dplyr)
data1$Gene_name <- ""
for(i in 1:length(data1$Gene_Id)){
  id <- data2 %>% filter(Gene_Id == as.character(data1$Gene_Id[i]))
  data1$Gene_name[i] = as.character(id$Gene_name[1])
}
write.csv(data1, file="data7.csv", row.names = F)

#ensembl
setwd("/home/dell/Documents/koushiki/project_new/ensembl/")

data1 <- read.csv("data6.csv", header=T)
head(data1)
data2 <- read.csv("gene_name.txt", header=F, sep = "\t")
head(data2)
colnames(data2) <- c("Gene_Id", "Gene_name")
library(dplyr)
data1$Gene_name <- ""
for(i in 1:length(data1$Gene_Id)){
  id <- data2 %>% filter(Gene_Id == as.character(data1$Gene_Id[i]))
  data1$Gene_name[i] = as.character(id$Gene_name[1])
}
write.csv(data1, file="data7.csv", row.names = F)
