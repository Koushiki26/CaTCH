setwd("/home/dell/Documents/koushiki/project/")

transcriptid <- read.table("transcript_id.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
exon_data <- read.table("exon_length1.txt", header = F)
colnames(exon_data) <- c("TranscriptID","Exon_Length")
head(exon_data)

library(dplyr)
exon_data$GeneID <- ''
for(i in 1:length(exon_data$TranscriptID)){
  trans <- transcriptid %>% filter(as.character(TranscriptID)==as.character(exon_data$TranscriptID[i]))
  exon_data$GeneID[i] <- as.character(trans$GeneID[1])
}
exon_data <- na.omit(exon_data)
write.csv(exon_data, file = "data_exon1.csv", col.names = T, row.names = F)

main_data <- read.csv("data5.csv", header = T)
#colnames(main_data) <- c("geneID","Gene_Length","Transcript","Exon")
head(main_data)
exon_data <- read.csv("data_exon1.csv", header = T)
head(exon_data)
exon_data$Exon_Length <- as.numeric(exon_data$Exon_Length)
summary(exon_data)

main_data$Exon_Length <- ''
for(i in 1:length(main_data$Gene_Id))
{
  exon <- exon_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$Exon_Length[i] = median(exon$Exon_Length)
}
write.csv(main_data, file = "data7.csv", col.names = T, row.names = F)
