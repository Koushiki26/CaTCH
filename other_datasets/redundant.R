setwd("/home/dell/Documents/koushiki/project_new/")

#to remove redundant genes, transcripts, and exons
data <- read.csv(file.choose(), header = F, sep = "\t")
colnames(data) <- c('seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute')
head(data)
library(dplyr)
x <- data %>% filter(feature=="gene")
head(x)
count(x)
y <- data %>% filter(feature=="gene") %>% distinct(start, end, .keep_all = T)
count(y)
a <- data %>% filter(feature=="exon")
count(a)
b <- data %>% filter(feature=="exon") %>% distinct(start, end, .keep_all = T)
head(b)
count(b)
c <- data %>% filter(feature=="transcript")
count(c)
d <- data %>% filter(feature=="transcript") %>% distinct(start, end, .keep_all=T)
head(d)
count(d)

#for GENCODE
write.table(y, file="gencode/gene.txt", sep="\t", quote = F, row.names =F, col.names = F )
write.table(d, file="gencode/transcript.txt", sep="\t", quote = F, row.names =F, col.names = F )
write.table(b, file="gencode/exon.txt", sep="\t", quote = F, row.names =F, col.names = F )

#for ensembl
write.table(y, file="ensembl/gene.txt", sep="\t", quote = F, row.names =F, col.names = F )
write.table(d, file="ensembl/transcript.txt", sep="\t", quote = F, row.names =F, col.names = F )
write.table(b, file="ensembl/exon.txt", sep="\t", quote = F, row.names =F, col.names = F )

#for refseq_T2T
write.table(y, file="refseq_T2T/gene.txt", sep="\t", quote = F, row.names =F, col.names = F )
write.table(d, file="refseq_T2T/transcript.txt", sep="\t", quote = F, row.names =F, col.names = F )
write.table(b, file="refseq_T2T/exon.txt", sep="\t", quote = F, row.names =F, col.names = F )