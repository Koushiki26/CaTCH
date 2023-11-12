#gencode
setwd("/home/dell/Documents/koushiki/project_new/gencode/")
library(dplyr)

cp_data <- read.csv("data_cp.txt", header = F, sep="\t")
colnames(cp_data) <- c("TranscriptId","GeneId","CP")
head(cp_data)

main_data <- read.csv("data4.csv", header = T)
head(main_data)
main_data$codingpotential <- ''
for(i in 1:length(main_data$Gene_Id)){
  cp <- cp_data %>% filter(as.character(GeneId) == as.character(main_data$Gene_Id[i]))
  main_data$codingpotential[i] = median(cp$CP)
}
write.csv(main_data, file = "data5.csv", col.names = T, row.names = F)

#ensembl
setwd("/home/dell/Documents/koushiki/project_new/ensembl/")
library(dplyr)

cp_data <- read.csv("data_cp.txt", header = F, sep="\t")
colnames(cp_data) <- c("TranscriptId","GeneId","CP")
head(cp_data)

main_data <- read.csv("data4.csv", header = T)
head(main_data)
main_data$codingpotential <- ''
for(i in 1:length(main_data$Gene_Id)){
  cp <- cp_data %>% filter(as.character(GeneId) == as.character(main_data$Gene_Id[i]))
  main_data$codingpotential[i] = median(cp$CP)
}
write.csv(main_data, file = "data5.csv", col.names = T, row.names = F)

#refseq_T2T
setwd("/home/dell/Documents/koushiki/project_new/refseq_T2T/")
library(dplyr)

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
cp_data <- read.table("data_cp_old.txt", header = F)
colnames(cp_data) <- c("TranscriptID","CP")
head(cp_data)

library(dplyr)
cp_data$GeneID <- ''
for(i in 1:length(cp_data$TranscriptID)){
  cp <- transcriptid %>% filter(as.character(TranscriptID)==as.character(cp_data$TranscriptID[i]))
  cp_data$GeneID[i] <- as.character(cp$GeneID[1])
}
cp_data <- na.omit(cp_data)
write.csv(cp_data, file = "data_cp.csv", col.names = T, row.names = F)

cp_data <- read.csv("data_cp.csv", header = T)
colnames(cp_data) <- c("TranscriptId","CP","GeneId")
head(cp_data)

main_data <- read.csv("data4.csv", header = T)
head(main_data)
main_data$codingpotential <- ''
for(i in 1:length(main_data$Gene_Id)){
  cp <- cp_data %>% filter(as.character(GeneId) == as.character(main_data$Gene_Id[i]))
  main_data$codingpotential[i] = median(cp$CP)
}
write.csv(main_data, file = "data5.csv", col.names = T, row.names = F)
