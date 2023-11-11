setwd("/home/dell/Documents/koushiki/project/")

data1 <- read.csv("data2.csv", header=F)
head(data1)
colnames(data1) <- c("Gene_Id", "Gene_length", "Num_of_Transcripts", "Num_of_Exons")
data2 <- read.csv("data3.csv", header=F)
head(data2)
colnames(data2) <- c("Gene_Id", "Gene_name", "Gene_type")
library(dplyr)
data1[,c("Gene_name", "Gene_type")] <- ""
for(i in 1:length(data1$Gene_Id)){
  id <- data2 %>% filter(Gene_Id == as.character(data1$Gene_Id[i]))
  data1$Gene_name[i] = as.character(id$Gene_name[1])
  data1$Gene_type[i] = as.character(id$Gene_type[1])
}
write.csv(data1, file="data4.csv", row.names = F)
