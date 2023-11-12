setwd("/home/dell/Documents/koushiki/project_new/refseq_T2T/")

data1 <- read.csv("gencode.v43.metadata.RefSeq", header = F, sep="\t")
data1 <- data1[,-3]
colnames(data1) <- c("Transcript1","Transcript2")
head(data1)
sss5_data <- read.table("sss5_1.txt", header = F, sep="\t")
colnames(sss5_data) <- c("TranscriptId", "strength")
head(sss5_data)
sss3_data <- read.table("sss3_1.txt", header = F, sep="\t")
colnames(sss3_data) <- c("TranscriptId", "strength")
head(sss3_data)

data1$strength5 <- ""
data1$strength3 <- ""
for(i in 1:nrow(data1)){
  s5 <- sss5_data %>% filter(as.character(TranscriptId)==as.character(data1$Transcript1[i]))
  data1$strength5[i] <- as.numeric(s5$strength[1])
  s3 <- sss3_data %>% filter(as.character(TranscriptId)==as.character(data1$Transcript1[i]))
  data1$strength3[i] <- as.numeric(s3$strength[1])
}
data1 <- na.omit(data1)
ss5 <- data1[,c(2,3)]
ss3 <- data1[,c(2,4)]
write.table(ss5, file = "sss5.txt", col.names = F, row.names = F, sep="\t")
write.table(ss3, file = "sss3.txt", col.names = F, row.names = F, sep="\t")
