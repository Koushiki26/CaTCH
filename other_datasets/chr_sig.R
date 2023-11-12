#gencode
setwd("/home/dell/Documents/koushiki/project_new/gencode/")

chr_sig1 <- read.table("../H3K36me3_gene.txt", header=F)
colnames(chr_sig1) <- c("peakscore", "geneid","genename")
head(chr_sig1)
chr_sig2 <- read.table("../H3K9me3_gene.txt", header=F)
colnames(chr_sig2) <- c("peakscore", "geneid","genename")
head(chr_sig2)
chr_sig3 <- read.table("../H3K4me3_gene.txt", header=F)
colnames(chr_sig3) <- c("peakscore", "geneid","genename")
head(chr_sig3)

main_data <- read.csv("data7.csv", header = T)
head(main_data)
main_data$H3K36me3 <- ''
main_data$H3K9me3 <- ''
main_data$H3K4me3 <- ''
for(i in 1:nrow(main_data)){
  test = grepl(as.character(main_data$Gene_name[i]), as.character(chr_sig1$genename[1:nrow(chr_sig1)]))
  if(sum(test)>=1){
    main_data$H3K36me3[i] = 1
  }
  else{
    main_data$H3K36me3[i] = 0
  }
}
for(i in 1:nrow(main_data)){
  test = grepl(as.character(main_data$Gene_name[i]), as.character(chr_sig2$genename[1:nrow(chr_sig2)]))
  if(sum(test)>=1){
    main_data$H3K9me3[i] = 1
  }
  else{
    main_data$H3K9me3[i] = 0
  }
}
for(i in 1:nrow(main_data)){
  test = grepl(as.character(main_data$Gene_name[i]), as.character(chr_sig3$genename[1:nrow(chr_sig3)]))
  if(sum(test)>=1){
    main_data$H3K4me3[i] = 1
  }
  else{
    main_data$H3K4me3[i] = 0
  }
}
main_data$chr_sig <- ''
for(i in 1:nrow(main_data)){
  main_data$chr_sig[i] <- sum(as.numeric(main_data[i,13:15]))
}
write.csv(main_data, file = "data8.csv", col.names = T, row.names = F)

#ensembl
setwd("/home/dell/Documents/koushiki/project_new/ensembl/")

chr_sig1 <- read.table("../H3K36me3_gene.txt", header=F)
colnames(chr_sig1) <- c("peakscore", "geneid","genename")
head(chr_sig1)
chr_sig2 <- read.table("../H3K9me3_gene.txt", header=F)
colnames(chr_sig2) <- c("peakscore", "geneid","genename")
head(chr_sig2)
chr_sig3 <- read.table("../H3K4me3_gene.txt", header=F)
colnames(chr_sig3) <- c("peakscore", "geneid","genename")
head(chr_sig3)

main_data <- read.csv("data7.csv", header = T)
head(main_data)
main_data <- na.omit(main_data)
main_data$H3K36me3 <- ''
main_data$H3K9me3 <- ''
main_data$H3K4me3 <- ''
for(i in 1:nrow(main_data)){
  test = grepl(as.character(main_data$Gene_name[i]), as.character(chr_sig1$genename[1:nrow(chr_sig1)]))
  if(sum(test)>=1){
    main_data$H3K36me3[i] = 1
  }
  else{
    main_data$H3K36me3[i] = 0
  }
}
for(i in 1:nrow(main_data)){
  test = grepl(as.character(main_data$Gene_name[i]), as.character(chr_sig2$genename[1:nrow(chr_sig2)]))
  if(sum(test)>=1){
    main_data$H3K9me3[i] = 1
  }
  else{
    main_data$H3K9me3[i] = 0
  }
}
for(i in 1:nrow(main_data)){
  test = grepl(as.character(main_data$Gene_name[i]), as.character(chr_sig3$genename[1:nrow(chr_sig3)]))
  if(sum(test)>=1){
    main_data$H3K4me3[i] = 1
  }
  else{
    main_data$H3K4me3[i] = 0
  }
}
main_data$chr_sig <- ''
for(i in 1:nrow(main_data)){
  main_data$chr_sig[i] <- sum(as.numeric(main_data[i,13:15]))
}
write.csv(main_data, file = "data8.csv", col.names = T, row.names = F)

#refseq_T2T
setwd("/home/dell/Documents/koushiki/project_new/refseq_T2T/")

chr_sig1 <- read.table("../H3K36me3_gene.txt", header=F)
colnames(chr_sig1) <- c("peakscore", "geneid","genename")
head(chr_sig1)
chr_sig2 <- read.table("../H3K9me3_gene.txt", header=F)
colnames(chr_sig2) <- c("peakscore", "geneid","genename")
head(chr_sig2)
chr_sig3 <- read.table("../H3K4me3_gene.txt", header=F)
colnames(chr_sig3) <- c("peakscore", "geneid","genename")
head(chr_sig3)

main_data <- read.csv("data6.csv", header = T)
head(main_data)
main_data$H3K36me3 <- ''
main_data$H3K9me3 <- ''
main_data$H3K4me3 <- ''
for(i in 1:nrow(main_data)){
  test = grepl(as.character(main_data$Gene_Id[i]), as.character(chr_sig1$genename[1:nrow(chr_sig1)]))
  if(sum(test)>=1){
    main_data$H3K36me3[i] = 1
  }
  else{
    main_data$H3K36me3[i] = 0
  }
}
for(i in 1:nrow(main_data)){
  test = grepl(as.character(main_data$Gene_Id[i]), as.character(chr_sig2$genename[1:nrow(chr_sig2)]))
  if(sum(test)>=1){
    main_data$H3K9me3[i] = 1
  }
  else{
    main_data$H3K9me3[i] = 0
  }
}
for(i in 1:nrow(main_data)){
  test = grepl(as.character(main_data$Gene_Id[i]), as.character(chr_sig3$genename[1:nrow(chr_sig3)]))
  if(sum(test)>=1){
    main_data$H3K4me3[i] = 1
  }
  else{
    main_data$H3K4me3[i] = 0
  }
}
main_data$chr_sig <- ''
for(i in 1:nrow(main_data)){
  main_data$chr_sig[i] <- sum(as.numeric(main_data[i,12:14]))
}
write.csv(main_data, file = "data7.csv", col.names = T, row.names = F)