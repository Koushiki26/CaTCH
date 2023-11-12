#gencode
setwd("/home/dell/Documents/koushiki/project_new/gencode/")

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
ss_data <- read.table("splicesite.txt", header = F)
colnames(ss_data) <- c("TranscriptID","ss5","ss3")
head(ss_data)

library(dplyr)
ss_data$GeneID <- ''
for(i in 1:length(ss_data$TranscriptID)){
  ss <- transcriptid %>% filter(as.character(TranscriptID)==as.character(ss_data$TranscriptID[i]))
  ss_data$GeneID[i] <- as.character(ss$GeneID[1])
}
ss_data <- na.omit(ss_data)
write.csv(ss_data, file = "data_splicesite.csv", col.names = T, row.names = F)

main_data <- read.csv("data3.csv", header = T)
head(main_data)
ss_data <- read.csv("data_splicesite.csv", header = T)
head(ss_data)
ss_data$ss5 <- as.factor(ss_data$ss5)
ss_data$ss3 <- as.factor(ss_data$ss3)
summary(ss_data)

main_data$ss5 <- ''
main_data$ss3 <- ''
for(i in 1:length(main_data$Gene_Id))
{
  splice <- ss_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$ss5[i] = list(unique(as.character(splice$ss5)))
  main_data$ss3[i] = list(unique(as.character(splice$ss3)))
}
main_data <- na.omit(main_data)

x <- c("GT","GC","AT")
comb1 <- combn(x,1)
comb2 <- combn(x,2)
comb3 <- combn(x,3)

main_data$GT <- ''
for(i in 1:nrow(main_data)){
  test1 <- grepl("GT", main_data$ss5[i])
  if(test1==T){
    main_data$GT[i] = 1
  }
  else{
    main_data$GT[i] = 0
  }
}
main_data$GC <- ''
for(i in 1:nrow(main_data)){
  test2 <- grepl("GC", main_data$ss5[i])
  if(test2==T){
    main_data$GC[i] = 2
  }
  else{
    main_data$GC[i] = 0
  }
}
main_data$AT <- ''
for(i in 1:nrow(main_data)){
  test3 <- grepl("AT", main_data$ss5[i])
  if(test3==T){
    main_data$AT[i] = 3
  }
  else{
    main_data$AT[i] = 0
  }
}
main_data$GTGC <- ''
for(i in 1:nrow(main_data)){
  test4 <- grepl("GT", main_data$ss5[i]) & grepl("GC", main_data$ss5[i])
  if(test4==T){
    main_data$GTGC[i] = replace(main_data$splice5[i],1,4)
  }
  else{
    main_data$GTGC[i] = 0
  }
}
main_data$GTAT <- ''
for(i in 1:nrow(main_data)){
  test5 <- grepl("GT", main_data$ss5[i]) & grepl("AT", main_data$ss5[i])
  if(test5==T){
    main_data$GTAT[i] = 5
  }
  else{
    main_data$GTAT[i] = 0
  }
}
main_data$GCAT <- ''
for(i in 1:nrow(main_data)){
  test6 <- grepl("GC", main_data$ss5[i]) & grepl("AT", main_data$ss5[i])
  if(test6==T){
    main_data$GCAT[i] = 6
  }
  else{
    main_data$GCAT[i] = 0
  }
}
main_data$GTGCAT <- ''
for(i in 1:nrow(main_data)){
  test7 <- grepl("GT", main_data$ss5[i]) & grepl("GC", main_data$ss5[i]) & grepl("AT", main_data$ss5[i])
  if(test7==T){
    main_data$GTGCAT[i] = 7
  }
  else{
    main_data$GTGCAT[i] = 0
  }
}
main_data$splice5 <- ''
for(i in 1:nrow(main_data)){
  main_data$splice5[i] = max(main_data$GT[i],main_data$GC[i],main_data$AT[i],main_data$GTGC[i],
                             main_data$GTAT[i],main_data$GCAT[i],main_data$GTGCAT[i])
}

x <- c("AC","AG")
comb1 <- combn(x,1)
comb2 <- combn(x,2)

main_data$AC <- ''
for(i in 1:nrow(main_data)){
  test1 <- grepl("AC", main_data$ss5[i])
  if(test1==T){
    main_data$AC[i] = 1
  }
  else{
    main_data$AC[i] = 0
  }
}
main_data$AG <- ''
for(i in 1:nrow(main_data)){
  test2 <- grepl("AG", main_data$ss5[i])
  if(test2==T){
    main_data$AG[i] = 2
  }
  else{
    main_data$AG[i] = 0
  }
}
main_data$ACAG <- ''
for(i in 1:nrow(main_data)){
  test3 <- grepl("AC", main_data$ss5[i]) & grepl("AG", main_data$ss5[i])
  if(test3==T){
    main_data$ACAG[i] = 3
  }
  else{
    main_data$ACAG[i] = 0
  }
}
main_data$splice3 <- ''
for(i in 1:nrow(main_data)){
  main_data$splice3[i] = max(main_data$AC[i],main_data$AG[i],main_data$ACAG[i])
}

names(main_data)
main_data1 <- main_data[,-c(8,9,10,11,12,13,14,15,16,18,19,20)]
write.csv(main_data1, file = "data4.csv", col.names = T, row.names = F)

#ensembl
setwd("/home/dell/Documents/koushiki/project_new/ensembl/")

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
ss_data <- read.table("splicesite.txt", header = F)
colnames(ss_data) <- c("TranscriptID","ss5","ss3")
head(ss_data)

library(dplyr)
ss_data$GeneID <- ''
for(i in 1:length(ss_data$TranscriptID)){
  ss <- transcriptid %>% filter(as.character(TranscriptID)==as.character(ss_data$TranscriptID[i]))
  ss_data$GeneID[i] <- as.character(ss$GeneID[1])
}
ss_data <- na.omit(ss_data)
write.csv(ss_data, file = "data_splicesite.csv", col.names = T, row.names = F)

main_data <- read.csv("data3.csv", header = T)
head(main_data)
ss_data <- read.csv("data_splicesite.csv", header = T)
head(ss_data)
ss_data$ss5 <- as.factor(ss_data$ss5)
ss_data$ss3 <- as.factor(ss_data$ss3)
summary(ss_data)

main_data$ss5 <- ''
main_data$ss3 <- ''
for(i in 1:length(main_data$Gene_Id))
{
  splice <- ss_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$ss5[i] = list(unique(as.character(splice$ss5)))
  main_data$ss3[i] = list(unique(as.character(splice$ss3)))
}
main_data <- na.omit(main_data)

x <- c("GT","GC","AT")
comb1 <- combn(x,1)
comb2 <- combn(x,2)
comb3 <- combn(x,3)

main_data$GT <- ''
for(i in 1:nrow(main_data)){
  test1 <- grepl("GT", main_data$ss5[i])
  if(test1==T){
    main_data$GT[i] = 1
  }
  else{
    main_data$GT[i] = 0
  }
}
main_data$GC <- ''
for(i in 1:nrow(main_data)){
  test2 <- grepl("GC", main_data$ss5[i])
  if(test2==T){
    main_data$GC[i] = 2
  }
  else{
    main_data$GC[i] = 0
  }
}
main_data$AT <- ''
for(i in 1:nrow(main_data)){
  test3 <- grepl("AT", main_data$ss5[i])
  if(test3==T){
    main_data$AT[i] = 3
  }
  else{
    main_data$AT[i] = 0
  }
}
main_data$GTGC <- ''
for(i in 1:nrow(main_data)){
  test4 <- grepl("GT", main_data$ss5[i]) & grepl("GC", main_data$ss5[i])
  if(test4==T){
    main_data$GTGC[i] = replace(main_data$splice5[i],1,4)
  }
  else{
    main_data$GTGC[i] = 0
  }
}
main_data$GTAT <- ''
for(i in 1:nrow(main_data)){
  test5 <- grepl("GT", main_data$ss5[i]) & grepl("AT", main_data$ss5[i])
  if(test5==T){
    main_data$GTAT[i] = 5
  }
  else{
    main_data$GTAT[i] = 0
  }
}
main_data$GCAT <- ''
for(i in 1:nrow(main_data)){
  test6 <- grepl("GC", main_data$ss5[i]) & grepl("AT", main_data$ss5[i])
  if(test6==T){
    main_data$GCAT[i] = 6
  }
  else{
    main_data$GCAT[i] = 0
  }
}
main_data$GTGCAT <- ''
for(i in 1:nrow(main_data)){
  test7 <- grepl("GT", main_data$ss5[i]) & grepl("GC", main_data$ss5[i]) & grepl("AT", main_data$ss5[i])
  if(test7==T){
    main_data$GTGCAT[i] = 7
  }
  else{
    main_data$GTGCAT[i] = 0
  }
}
main_data$splice5 <- ''
for(i in 1:nrow(main_data)){
  main_data$splice5[i] = max(main_data$GT[i],main_data$GC[i],main_data$AT[i],main_data$GTGC[i],
                             main_data$GTAT[i],main_data$GCAT[i],main_data$GTGCAT[i])
}

x <- c("AC","AG")
comb1 <- combn(x,1)
comb2 <- combn(x,2)

main_data$AC <- ''
for(i in 1:nrow(main_data)){
  test1 <- grepl("AC", main_data$ss5[i])
  if(test1==T){
    main_data$AC[i] = 1
  }
  else{
    main_data$AC[i] = 0
  }
}
main_data$AG <- ''
for(i in 1:nrow(main_data)){
  test2 <- grepl("AG", main_data$ss5[i])
  if(test2==T){
    main_data$AG[i] = 2
  }
  else{
    main_data$AG[i] = 0
  }
}
main_data$ACAG <- ''
for(i in 1:nrow(main_data)){
  test3 <- grepl("AC", main_data$ss5[i]) & grepl("AG", main_data$ss5[i])
  if(test3==T){
    main_data$ACAG[i] = 3
  }
  else{
    main_data$ACAG[i] = 0
  }
}
main_data$splice3 <- ''
for(i in 1:nrow(main_data)){
  main_data$splice3[i] = max(main_data$AC[i],main_data$AG[i],main_data$ACAG[i])
}

names(main_data)
main_data1 <- main_data[,-c(8,9,10,11,12,13,14,15,16,18,19,20)]
write.csv(main_data1, file = "data4.csv", col.names = T, row.names = F)

#refseq_T2T
setwd("/home/dell/Documents/koushiki/project_new/refseq_T2T/")

transcriptid <- read.table("TranscriptId.txt", header = F)
colnames(transcriptid) <- c("GeneID","TranscriptID")
head(transcriptid)
ss_data <- read.table("splicesite.txt", header = F)
colnames(ss_data) <- c("TranscriptID","ss5","ss3")
head(ss_data)

library(dplyr)
ss_data$GeneID <- ''
for(i in 1:length(ss_data$TranscriptID)){
  ss <- transcriptid %>% filter(as.character(TranscriptID)==as.character(ss_data$TranscriptID[i]))
  ss_data$GeneID[i] <- as.character(ss$GeneID[1])
}
ss_data <- na.omit(ss_data)
write.csv(ss_data, file = "data_splicesite.csv", col.names = T, row.names = F)

main_data <- read.csv("data3.csv", header = T)
head(main_data)
ss_data <- read.csv("data_splicesite.csv", header = T)
head(ss_data)
ss_data$ss5 <- as.factor(ss_data$ss5)
ss_data$ss3 <- as.factor(ss_data$ss3)
summary(ss_data)

main_data$ss5 <- ''
main_data$ss3 <- ''
for(i in 1:length(main_data$Gene_Id))
{
  splice <- ss_data %>% filter(as.character(GeneID) == as.character(main_data$Gene_Id[i]))
  main_data$ss5[i] = list(unique(as.character(splice$ss5)))
  main_data$ss3[i] = list(unique(as.character(splice$ss3)))
}
main_data <- na.omit(main_data)

x <- c("GT","GC","AT")
comb1 <- combn(x,1)
comb2 <- combn(x,2)
comb3 <- combn(x,3)

main_data$GT <- ''
for(i in 1:nrow(main_data)){
  test1 <- grepl("GT", main_data$ss5[i])
  if(test1==T){
    main_data$GT[i] = 1
  }
  else{
    main_data$GT[i] = 0
  }
}
main_data$GC <- ''
for(i in 1:nrow(main_data)){
  test2 <- grepl("GC", main_data$ss5[i])
  if(test2==T){
    main_data$GC[i] = 2
  }
  else{
    main_data$GC[i] = 0
  }
}
main_data$AT <- ''
for(i in 1:nrow(main_data)){
  test3 <- grepl("AT", main_data$ss5[i])
  if(test3==T){
    main_data$AT[i] = 3
  }
  else{
    main_data$AT[i] = 0
  }
}
main_data$GTGC <- ''
for(i in 1:nrow(main_data)){
  test4 <- grepl("GT", main_data$ss5[i]) & grepl("GC", main_data$ss5[i])
  if(test4==T){
    main_data$GTGC[i] = replace(main_data$splice5[i],1,4)
  }
  else{
    main_data$GTGC[i] = 0
  }
}
main_data$GTAT <- ''
for(i in 1:nrow(main_data)){
  test5 <- grepl("GT", main_data$ss5[i]) & grepl("AT", main_data$ss5[i])
  if(test5==T){
    main_data$GTAT[i] = 5
  }
  else{
    main_data$GTAT[i] = 0
  }
}
main_data$GCAT <- ''
for(i in 1:nrow(main_data)){
  test6 <- grepl("GC", main_data$ss5[i]) & grepl("AT", main_data$ss5[i])
  if(test6==T){
    main_data$GCAT[i] = 6
  }
  else{
    main_data$GCAT[i] = 0
  }
}
main_data$GTGCAT <- ''
for(i in 1:nrow(main_data)){
  test7 <- grepl("GT", main_data$ss5[i]) & grepl("GC", main_data$ss5[i]) & grepl("AT", main_data$ss5[i])
  if(test7==T){
    main_data$GTGCAT[i] = 7
  }
  else{
    main_data$GTGCAT[i] = 0
  }
}
main_data$splice5 <- ''
for(i in 1:nrow(main_data)){
  main_data$splice5[i] = max(main_data$GT[i],main_data$GC[i],main_data$AT[i],main_data$GTGC[i],
                             main_data$GTAT[i],main_data$GCAT[i],main_data$GTGCAT[i])
}

x <- c("AC","AG")
comb1 <- combn(x,1)
comb2 <- combn(x,2)

main_data$AC <- ''
for(i in 1:nrow(main_data)){
  test1 <- grepl("AC", main_data$ss5[i])
  if(test1==T){
    main_data$AC[i] = 1
  }
  else{
    main_data$AC[i] = 0
  }
}
main_data$AG <- ''
for(i in 1:nrow(main_data)){
  test2 <- grepl("AG", main_data$ss5[i])
  if(test2==T){
    main_data$AG[i] = 2
  }
  else{
    main_data$AG[i] = 0
  }
}
main_data$ACAG <- ''
for(i in 1:nrow(main_data)){
  test3 <- grepl("AC", main_data$ss5[i]) & grepl("AG", main_data$ss5[i])
  if(test3==T){
    main_data$ACAG[i] = 3
  }
  else{
    main_data$ACAG[i] = 0
  }
}
main_data$splice3 <- ''
for(i in 1:nrow(main_data)){
  main_data$splice3[i] = max(main_data$AC[i],main_data$AG[i],main_data$ACAG[i])
}

names(main_data)
main_data1 <- main_data[,-c(8,9,10,11,12,13,14,15,16,18,19,20)]
write.csv(main_data1, file = "data4.csv", col.names = T, row.names = F)
