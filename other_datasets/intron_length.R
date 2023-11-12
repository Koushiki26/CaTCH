#gencode
setwd("/home/dell/Documents/koushiki/project_new/gencode/")

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
intron_data <- read.table("intron_length.txt", header = F)
colnames(intron_data) <- c("TranscriptID","Intron_Length")
head(intron_data)

library(dplyr)
intron_data$GeneID <- ''
for(i in 1:length(intron_data$TranscriptID)){
  trans <- transcriptid %>% filter(as.character(TranscriptID)==as.character(intron_data$TranscriptID[i]))
  intron_data$GeneID[i] <- as.character(trans$GeneID[1])
}
intron_data <- na.omit(intron_data)
write.csv(intron_data, file = "data_intron.csv", col.names = T, row.names = F)

main_data <- read.csv("data2.csv", header = T)
head(main_data)
intron_data <- read.csv("data_intron.csv", header = T)
head(intron_data)
intron_data$Intron_Length <- as.numeric(intron_data$Intron_Length)
summary(intron_data)

main_data$Intron_Length <- ''
for(i in 1:length(main_data$Gene_Id))
{
  intron <- intron_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$Intron_Length[i] = median(intron$Intron_Length)
}
write.csv(main_data, file = "data3.csv", col.names = T, row.names = F)

#ensembl
setwd("/home/dell/Documents/koushiki/project_new/ensembl/")

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
intron_data <- read.table("intron_length.txt", header = F)
colnames(intron_data) <- c("TranscriptID","Intron_Length")
head(intron_data)

library(dplyr)
intron_data$GeneID <- ''
for(i in 1:length(intron_data$TranscriptID)){
  trans <- transcriptid %>% filter(as.character(TranscriptID)==as.character(intron_data$TranscriptID[i]))
  intron_data$GeneID[i] <- as.character(trans$GeneID[1])
}
intron_data <- na.omit(intron_data)
write.csv(intron_data, file = "data_intron.csv", col.names = T, row.names = F)

main_data <- read.csv("data2.csv", header = T)
head(main_data)
intron_data <- read.csv("data_intron.csv", header = T)
head(intron_data)
intron_data$Intron_Length <- as.numeric(intron_data$Intron_Length)
summary(intron_data)

main_data$Intron_Length <- ''
for(i in 1:length(main_data$Gene_Id))
{
  intron <- intron_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$Intron_Length[i] = median(intron$Intron_Length)
}
write.csv(main_data, file = "data3.csv", col.names = T, row.names = F)

#refseq_T2T
setwd("/home/dell/Documents/koushiki/project_new/refseq_T2T/")

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
intron_data <- read.table("intron_length.txt", header = F)
colnames(intron_data) <- c("TranscriptID","Intron_Length")
head(intron_data)

library(dplyr)
intron_data$GeneID <- ''
for(i in 1:length(intron_data$TranscriptID)){
  trans <- transcriptid %>% filter(as.character(TranscriptID)==as.character(intron_data$TranscriptID[i]))
  intron_data$GeneID[i] <- as.character(trans$GeneID[1])
}
intron_data <- na.omit(intron_data)
write.csv(intron_data, file = "data_intron.csv", col.names = T, row.names = F)

main_data <- read.csv("data2.csv", header = T)
head(main_data)
intron_data <- read.csv("data_intron.csv", header = T)
head(intron_data)
intron_data$Intron_Length <- as.numeric(intron_data$Intron_Length)
summary(intron_data)

main_data$Intron_Length <- ''
for(i in 1:length(main_data$Gene_Id))
{
  intron <- intron_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$Intron_Length[i] = median(intron$Intron_Length)
}
write.csv(main_data, file = "data3.csv", col.names = T, row.names = F)
