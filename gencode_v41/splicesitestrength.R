setwd("/home/dell/Documents/koushiki/project/")

transcriptid <- read.table("transcript_id.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
sss5 <- read.table("/home/dell/Documents/koushiki/project/splicesitestrength/splicesite5strength.txt", header = F, sep="\t")
colnames(sss5) <- c("TranscriptId","strength5")
head(sss5)
sss3 <- read.table("/home/dell/Documents/koushiki/project/splicesitestrength/splicesite3strength.txt", header = F, sep="\t")
colnames(sss3) <- c("TranscriptId","strength3")
head(sss3)

library(dplyr)
sss5$GeneId <- ''
for(i in 1:length(sss5$TranscriptId)){
  strength <- transcriptid %>% filter(as.character(TranscriptID)==as.character(sss5$TranscriptId[i]))
  sss5$GeneId[i] <- as.character(strength$GeneID[1])
}
sss5 <- na.omit(sss5)
write.csv(sss5, file = "data_5sss.csv", col.names = T, row.names = F)

sss3$GeneId <- ''
for(i in 1:length(sss3$TranscriptId)){
  strength <- transcriptid %>% filter(as.character(TranscriptID)==as.character(sss3$TranscriptId[i]))
  sss3$GeneId[i] <- as.character(strength$GeneID[1])
}
sss3 <- na.omit(sss3)
write.csv(sss3, file = "data_3sss.csv", col.names = T, row.names = F)

main_data <- read.csv("data8.csv", header = T)
#colnames(main_data) <- c("geneID","Gene_Length","Transcript","Exon")
head(main_data)
sss5_data <- read.csv("data_5sss.csv", header = T)
head(sss5_data)
sss3_data <- read.csv("data_3sss.csv", header = T)
head(sss3_data)
main_data$ssstrength5 <- ''
main_data$ssstrength3 <- ''
for(i in 1:length(main_data$Gene_Id))
{
  splice5 <- sss5_data %>% filter(as.character(GeneId) == as.character(main_data$Gene_Id[i]))
  main_data$ssstrength5[i] = median(splice5$strength5)
  splice3 <- sss3_data %>% filter(as.character(GeneId) == as.character(main_data$Gene_Id[i]))
  main_data$ssstrength3[i] = median(splice3$strength3)
}
write.csv(main_data, file = "data9.csv", col.names = T, row.names = F)
