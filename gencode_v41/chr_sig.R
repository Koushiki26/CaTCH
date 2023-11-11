setwd("/home/dell/Documents/koushiki/project/")

chr_sig1 <- read.table("chrm_sig/new/H3K36me3_outgene.txt", header=F)
colnames(chr_sig1) <- c("peakscore", "geneid","genename")
head(chr_sig1)
chr_sig2 <- read.table("chrm_sig/new/H3K9me3_outgene.txt", header=F)
colnames(chr_sig2) <- c("peakscore", "geneid","genename")
head(chr_sig2)
chr_sig3 <- read.table("chrm_sig/new/H3K4me3_outgene.txt", header=F)
colnames(chr_sig3) <- c("peakscore", "geneid","genename")
head(chr_sig3)

main_data <- read.csv("data14.csv", header = T)
head(main_data)
main_data <- main_data[,-c(17,18,19)]
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
main_data$chr_sig_new <- ''
for(i in 1:nrow(main_data)){
  main_data$chr_sig_new[i] <- sum(as.numeric(main_data[i,18:20]))
}
write.csv(main_data, file = "data15.csv", col.names = T, row.names = F)

for(i in 1:nrow(main_data)){
  if(all(main_data$chr_sig[i]==main_data$chr_sig_new[i])!=T){
    print(i)
  }
}