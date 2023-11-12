#gencode
setwd("/home/dell/Documents/koushiki/project_new/gencode/")

data1 <- read.csv("main_data.csv", header=F)
head(data1)
colnames(data1) <- c("Gene_Id", "Gene_length", "Num_of_Transcripts", "Num_of_Exons")
data2 <- read.csv("gene_type.txt", header=F, sep = "\t")
head(data2)
colnames(data2) <- c("Gene_Id", "Gene_type")
library(dplyr)
data1$Gene_type <- ""
for(i in 1:length(data1$Gene_Id)){
  id <- data2 %>% filter(Gene_Id == as.character(data1$Gene_Id[i]))
  data1$Gene_type[i] = as.character(id$Gene_type[1])
}
write.csv(data1, file="data1.csv", row.names = F)

#ensembl
setwd("/home/dell/Documents/koushiki/project_new/ensembl/")

data1 <- read.csv("main_data.csv", header=F)
head(data1)
colnames(data1) <- c("Gene_Id", "Gene_length", "Num_of_Transcripts", "Num_of_Exons")
data2 <- read.csv("gene_type.txt", header=F, sep = "\t")
head(data2)
colnames(data2) <- c("Gene_Id", "Gene_type")
library(dplyr)
data1$Gene_type <- ""
for(i in 1:length(data1$Gene_Id)){
  id <- data2 %>% filter(Gene_Id == as.character(data1$Gene_Id[i]))
  data1$Gene_type[i] = as.character(id$Gene_type[1])
}
write.csv(data1, file="data1.csv", row.names = F)

#refseq_T2T
setwd("/home/dell/Documents/koushiki/project_new/refseq_T2T/")

data1 <- read.csv("main_data.csv", header=F)
head(data1)
colnames(data1) <- c("Gene_Id", "Gene_length", "Num_of_Transcripts", "Num_of_Exons")
data2 <- read.csv("gene_type.txt", header=F, sep = "\t")
head(data2)
colnames(data2) <- c("Gene_Id", "Gene_type")
library(dplyr)
data1$Gene_type <- ""
for(i in 1:length(data1$Gene_Id)){
  id <- data2 %>% filter(Gene_Id == as.character(data1$Gene_Id[i]))
  data1$Gene_type[i] = as.character(id$Gene_type[1])
}
write.csv(data1, file="data1.csv", row.names = F)
