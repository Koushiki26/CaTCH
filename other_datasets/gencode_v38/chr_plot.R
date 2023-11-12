setwd("/home/dell/Documents/koushiki/project_new/gencode/")
library(dplyr)

#chromatin_signature
main_data <- read.csv("data9.csv", header = T)
head(main_data)
main_data$TC <- main_data$Num_of_Transcripts/main_data$Num_of_Exons
main_data <- na.omit(main_data)
main_data <- main_data[,c(1,5,13,14,15,19)]
main_data$chr4_9 <- ''
for(i in 1:nrow(main_data)){
  main_data$chr4_9[i] <- sum(as.numeric(main_data[i,4:5]))
}
main_data$chr36_9 <- ''
for(i in 1:nrow(main_data)){
  main_data$chr36_9[i] <- sum(as.numeric(main_data[i,3:4]))
}
main_data$chr4_36 <- ''
for(i in 1:nrow(main_data)){
  main_data$chr4_36[i] <- sum(as.numeric(main_data[i,c(3,5)]))
}
main_data$chr4_9_36 <- ''
for(i in 1:nrow(main_data)){
  main_data$chr4_9_36[i] <- sum(as.numeric(main_data[i,3:5]))
}

lncrna_data <- main_data %>% filter(Gene_type=="lncRNA")
lncrna_data_low <- lncrna_data %>% filter(TC<=0.3333333)
lncrna_data_high <- lncrna_data %>% filter(TC>0.3333333)
mrna_data <- main_data %>% filter(Gene_type=="protein_coding")
mrna_data_low <- mrna_data %>% filter(TC<=0.3333333)
mrna_data_high <- mrna_data %>% filter(TC>0.3333333)

k36_1 <- lncrna_data_low %>% filter(H3K36me3==1)
a1 <- as.numeric(count(k36_1)/count(lncrna_data_low))
k36_2 <- lncrna_data_high %>% filter(H3K36me3==1)
b1 <- as.numeric(count(k36_2)/count(lncrna_data_high))
k36_3 <- mrna_data_low %>% filter(H3K36me3==1)
c1 <- as.numeric(count(k36_3)/count(mrna_data_low))
k36_4 <- mrna_data_high %>% filter(H3K36me3==1)
d1 <- as.numeric(count(k36_4)/count(mrna_data_high))
k9_1 <- lncrna_data_low %>% filter(H3K9me3==1)
a2 <- as.numeric(count(k9_1)/count(lncrna_data_low))
k9_2 <- lncrna_data_high %>% filter(H3K9me3==1)
b2 <- as.numeric(count(k9_2)/count(lncrna_data_high))
k9_3 <- mrna_data_low %>% filter(H3K9me3==1)
c2 <- as.numeric(count(k9_3)/count(mrna_data_low))
k9_4 <- mrna_data_high %>% filter(H3K9me3==1)
d2 <- as.numeric(count(k9_4)/count(mrna_data_high))
k4_1 <- lncrna_data_low %>% filter(H3K4me3==1)
a3 <- as.numeric(count(k4_1)/count(lncrna_data_low))
k4_2 <- lncrna_data_high %>% filter(H3K4me3==1)
b3 <- as.numeric(count(k4_2)/count(lncrna_data_high))
k4_3 <- mrna_data_low %>% filter(H3K4me3==1)
c3 <- as.numeric(count(k4_3)/count(mrna_data_low))
k4_4 <- mrna_data_high %>% filter(H3K4me3==1)
d3 <- as.numeric(count(k4_4)/count(mrna_data_high))
k36_4_1 <- lncrna_data_low %>% filter(chr4_36==2)
a4 <- as.numeric(count(k36_4_1)/count(lncrna_data_low))
k36_4_2 <- lncrna_data_high %>% filter(chr4_36==2)
b4 <- as.numeric(count(k36_4_2)/count(lncrna_data_high))
k36_4_3 <- mrna_data_low %>% filter(chr4_36==2)
c4 <- as.numeric(count(k36_4_3)/count(mrna_data_low))
k36_4_4 <- mrna_data_high %>% filter(chr4_36==2)
d4 <- as.numeric(count(k36_4_4)/count(mrna_data_high))
k9_4_1 <- lncrna_data_low %>% filter(chr4_9==2)
a5 <- as.numeric(count(k9_4_1)/count(lncrna_data_low))
k9_4_2 <- lncrna_data_high %>% filter(chr4_9==2)
b5 <- as.numeric(count(k9_4_2)/count(lncrna_data_high))
k9_4_3 <- mrna_data_low %>% filter(chr4_9==2)
c5 <- as.numeric(count(k9_4_3)/count(mrna_data_low))
k9_4_4 <- mrna_data_high %>% filter(chr4_9==2)
d5 <- as.numeric(count(k9_4_4)/count(mrna_data_high))
k9_36_1 <- lncrna_data_low %>% filter(chr36_9==2)
a6 <- as.numeric(count(k9_36_1)/count(lncrna_data_low))
k9_36_2 <- lncrna_data_high %>% filter(chr36_9==2)
b6 <- as.numeric(count(k9_36_2)/count(lncrna_data_high))
k9_36_3 <- mrna_data_low %>% filter(chr36_9==2)
c6 <- as.numeric(count(k9_36_3)/count(mrna_data_low))
k9_36_4 <- mrna_data_high %>% filter(chr36_9==2)
d6 <- as.numeric(count(k9_36_4)/count(mrna_data_high))
k9_36_4_1 <- lncrna_data_low %>% filter(chr4_9_36==3)
a7 <- as.numeric(count(k9_36_4_1)/count(lncrna_data_low))
k9_36_4_2 <- lncrna_data_high %>% filter(chr4_9_36==3)
b7 <- as.numeric(count(k9_36_4_2)/count(lncrna_data_high))
k9_36_4_3 <- mrna_data_low %>% filter(chr4_9_36==3)
c7 <- as.numeric(count(k9_36_4_3)/count(mrna_data_low))
k9_36_4_4 <- mrna_data_high %>% filter(chr4_9_36==3)
d7 <- as.numeric(count(k9_36_4_4)/count(mrna_data_high))

chrsig <- c(rep("H3K36me3",4),rep("H3K9me3",4),rep("H3K4me3",4),rep("H3K36me3 & H3K9me3",4),rep("H3K36me3 & H3K4me3",4),rep("H3K9me3 & H3K4me3",4),rep("all",4))
Type_of_gene <- rep(c("lncRNA (less TPE)","lncRNA (more TPE)","mRNA (less TPE)","mRNA (more TPE)"),7)
Frequency_of_genes <- c(a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,a6,b6,c6,d6,a4,b4,c4,d4,a5,b5,c5,d5,a7,b7,c7,d7)
df <- data.frame(chrsig,Type_of_gene,Frequency_of_genes)
df
level_order <- c("H3K36me3","H3K9me3","H3K4me3","H3K36me3 & H3K9me3","H3K36me3 & H3K4me3","H3K9me3 & H3K4me3","all")
tiff("/home/dell/Documents/koushiki/project_new/gencode/chr.tiff",width = 2400, height = 1200, res = 300)
p <- ggplot(df, aes(factor(chrsig, level=level_order), Frequency_of_genes, fill=Type_of_gene)) +
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = c("#31A2AC","#3A5199","#EB8A3E","#CB0000"))
p + xlab("") + ylab("Fraction of genes")+ 
  theme(axis.title = element_text(size = 19, color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black", face="bold")) + theme(legend.position = "none") +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
dev.off()