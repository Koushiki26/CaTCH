#gencode
setwd("/home/dell/Documents/koushiki/project_new/gencode/")

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
branch_data <- read.table("branch_count_length.txt", header = F)
colnames(branch_data) <- c("TranscriptID","count","min_length")
head(branch_data)

library(dplyr)
branch_data$GeneID <- ''
for(i in 1:length(branch_data$TranscriptID)){
  branch <- transcriptid %>% filter(as.character(TranscriptID)==as.character(branch_data$TranscriptID[i]))
  branch_data$GeneID[i] <- as.character(branch$GeneID[1])
}
branch_data <- na.omit(branch_data)
write.csv(branch_data, file = "data_branchpoint.csv", col.names = T, row.names = F)

branch_data <- read.csv("data_branchpoint.csv", header=T)
main_data <- read.csv("data5.csv", header = T)
head(main_data)
main_data$branchpt_len <- ''
for(i in 1:length(main_data$Gene_Id)){
  branch <- branch_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$branchpt_len[i] = median(branch$min_length, na.rm = FALSE)
}
write.csv(main_data, file = "data6.csv", col.names = T, row.names = F)

#ensembl
setwd("/home/dell/Documents/koushiki/project_new/ensembl/")

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
branch_data <- read.table("branch_count_length.txt", header = F)
colnames(branch_data) <- c("TranscriptID","count","min_length")
head(branch_data)

library(dplyr)
branch_data$GeneID <- ''
for(i in 1:length(branch_data$TranscriptID)){
  branch <- transcriptid %>% filter(as.character(TranscriptID)==as.character(branch_data$TranscriptID[i]))
  branch_data$GeneID[i] <- as.character(branch$GeneID[1])
}
branch_data <- na.omit(branch_data)
write.csv(branch_data, file = "data_branchpoint.csv", col.names = T, row.names = F)

branch_data <- read.csv("data_branchpoint.csv", header=T)
main_data <- read.csv("data5.csv", header = T)
head(main_data)
main_data$branchpt_len <- ''
for(i in 1:length(main_data$Gene_Id)){
  branch <- branch_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$branchpt_len[i] = median(branch$min_length, na.rm = FALSE)
}
write.csv(main_data, file = "data6.csv", col.names = T, row.names = F)

#refseq_T2T
setwd("/home/dell/Documents/koushiki/project_new/refseq_T2T/")

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
branch_data <- read.table("branch_count_length.txt", header = F)
colnames(branch_data) <- c("TranscriptID","count","min_length")
head(branch_data)

library(dplyr)
branch_data$GeneID <- ''
for(i in 1:length(branch_data$TranscriptID)){
  branch <- transcriptid %>% filter(as.character(TranscriptID)==as.character(branch_data$TranscriptID[i]))
  branch_data$GeneID[i] <- as.character(branch$GeneID[1])
}
branch_data <- na.omit(branch_data)
write.csv(branch_data, file = "data_branchpoint.csv", col.names = T, row.names = F)

branch_data <- read.csv("data_branchpoint.csv", header=T)
main_data <- read.csv("data5.csv", header = T)
head(main_data)
main_data$branchpt_len <- ''
for(i in 1:length(main_data$Gene_Id)){
  branch <- branch_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$branchpt_len[i] = median(branch$min_length, na.rm = FALSE)
}
write.csv(main_data, file = "data6.csv", col.names = T, row.names = F)
