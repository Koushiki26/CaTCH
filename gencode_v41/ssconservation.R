setwd("/home/dell/Documents/koushiki/project/")

transcriptid <- read.table("transcript_id.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
coor_data <- read.table("coordinate.txt", header = F)
colnames(coor_data) <- c("TranscriptID","chr","ss5","ss3")
head(coor_data)

library(dplyr)
coor_data$GeneID <- ''
for(i in 1:length(coor_data$TranscriptID)){
  coor <- transcriptid %>% filter(as.character(TranscriptID)==as.character(coor_data$TranscriptID[i]))
  coor_data$GeneID[i] <- as.character(coor$GeneID[1])
}
coor_data <- na.omit(coor_data)
write.csv(coor_data, file = "data_coordinate.csv", col.names = T, row.names = F)
summary(coor_data)


coor_data$conservation5 <- ''
coor_data$conservation3 <- ''
chr15 <- read.table("/home/dell/Documents/koushiki/project/conservation/chr1_f5.txt", header=F)
colnames(chr15) <- c("position","score")
head(chr15)
chr13 <- read.table("/home/dell/Documents/koushiki/project/conservation/chr1_f3.txt", header=F)
colnames(chr13) <- c("position","score")
head(chr13)
chr25 <- read.table("/home/dell/Documents/koushiki/project/conservation/chr2_f5.txt", header=F)
colnames(chr25) <- c("position","score")
head(chr25)
chr23 <- read.table("/home/dell/Documents/koushiki/project/conservation/chr2_f3.txt", header=F)
colnames(chr23) <- c("position","score")
head(chr23)
coor_data1 <- coor_data %>% filter(chr=="chr1")
for(i in 1:length(coor_data1$TranscriptID)){
  cons5 <- chr15 %>% filter(position==coor_data1$ss5[i])
  coor_data1$conservation5[i] <- cons5$score[1]
  cons3 <- chr13 %>% filter(position==coor_data1$ss3[i])
  coor_data1$conservation3[i] <- cons3$score[1]
}
coor_data2 <- coor_data %>% filter(chr=="chr2")
for(i in 1:length(coor_data2$TranscriptID)){
  cons5 <- chr25 %>% filter(position==coor_data2$ss5[i])
  coor_data2$conservation5[i] <- cons5$score[1]
  cons3 <- chr23 %>% filter(position==coor_data2$ss3[i])
  coor_data2$conservation3[i] <- cons3$score[1]
}
coor_datas <- rbind(coor_data1,coor_data2)
main_data <- read.csv("data9.csv", header = T)
head(main_data)
main_data$ss5conservation <- ''
main_data$ss3conservation <- ''
for(i in 1:length(main_data$Gene_Id)){
  cons <- coor_datas %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$ss5conservation[i] = median(cons$conservation5, na.rm = FALSE)
  main_data$ss3conservation[i] = median(cons$conservation3, na.rm = FALSE)
}
write.csv(main_data, file = "data10.csv", col.names = T, row.names = F)
