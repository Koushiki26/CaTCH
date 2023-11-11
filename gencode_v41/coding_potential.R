setwd("/home/dell/Documents/koushiki/project/")

cp_data <- read.csv("data_codingpotential.csv", header = F)
colnames(cp_data) <- c("TranscriptId","GeneId","CP","label")
head(cp_data)

main_data <- read.csv("data10.csv", header = T)
head(main_data)
main_data$codingpotential <- ''
for(i in 1:length(main_data$Gene_Id)){
  cp <- cp_data %>% filter(as.character(GeneId) == as.character(main_data$Gene_Id[i]))
  main_data$codingpotential[i] = median(cp$CP)
}
write.csv(main_data, file = "data11.csv", col.names = T, row.names = F)
