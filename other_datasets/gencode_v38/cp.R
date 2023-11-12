setwd("/home/dell/Documents/koushiki/project_new/gencode/")
library(dplyr)

#coding potential
main_data <- read.csv("data9.csv", header = T)
head(main_data)
main_data$TC <- main_data$Num_of_Transcripts/main_data$Num_of_Exons
main_data <- na.omit(main_data)
main_data <- main_data[,c(1,5,10,11,19)]

lncrna_data <- main_data %>% filter(Gene_type=="lncRNA")
lncrna_data_low <- lncrna_data %>% filter(TC<=0.3333333)
lncrna_data_high <- lncrna_data %>% filter(TC>0.3333333)
mrna_data <- main_data %>% filter(Gene_type=="protein_coding")
mrna_data_low <- mrna_data %>% filter(TC<=0.3333333)
mrna_data_high <- mrna_data %>% filter(TC>0.3333333)

tiff("/home/dell/Documents/koushiki/project_new/gencode/cp.tiff",width = 1200, height = 1200, res = 400)
par(mar=c(1,5,1,0.5), cex.lab = 1.2, font.lab=2)
boxplot(lncrna_data_low$codingpotential, lncrna_data_high$codingpotential, col = c("#31A2AC","#3A5199"), xlab="",xaxt='n', ylab="Coding Potential", boxwex=0.9,
        at=c(1,2), outline=F, xlim=c(0.5,4.5), ylim=c(0,1), main="", cex.axis=1.2, font.axis=2)
boxplot(mrna_data_low$codingpotential, mrna_data_high$codingpotential, add= TRUE,col = c("#EB8A3E","#CB0000"), boxwex=0.9, at=c(3,4), outline=F, axes=F)
box(col="black", lwd=2)
dev.off()

lncrna_data_low$type <- "lncRNA (low TC)"
lncrna_data_high$type <- "lncRNA (high TC)"
mrna_data_low$type <- "mRNA (low TC)"
mrna_data_high$type <- "mRNA (high TC)"
total_data <- data.frame(rbind(lncrna_data_low,lncrna_data_high,mrna_data_low,mrna_data_high))

tiff("/home/dell/Documents/koushiki/project_new/gencode/cp1.tiff",width = 1200, height = 1200, res = 300)
p <- ggplot(total_data, aes(x = codingpotential, y = type, fill = type)) +
  geom_density_ridges() + theme_ridges() + theme(legend.position = "none") + 
  scale_fill_manual(values = c("#3A5199","#31A2AC","#CB0000","#EB8A3E"))
p + xlab("Coding Potential") + ylab("")+ xlim(c(-0.35,1.25)) +
  theme(plot.margin = margin(0.25, -0.5, 0.25, -2.75, "cm"),
        axis.title.x = element_text(size = 19, color = "black", face = "bold", hjust=0.6),
        axis.text = element_text(size = 16, color = "black", face="bold"),
        plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.6),
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
dev.off()

