setwd("/home/dell/Documents/koushiki/project_new/gencode/")
library(dplyr)

#branching point plot
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

tiff("/home/dell/Documents/koushiki/project_new/gencode/branchpt.tiff",width = 1200, height = 1200, res = 400)
par(mar=c(1,5,0.5,0.5), cex.lab = 0.8, font.lab=2)
boxplot(lncrna_data_low$branchpt_len, lncrna_data_high$branchpt_len, col = c("#31A2AC","#3A5199"), xlab="",xaxt='n', ylab="Distance b/w branch point and 3'SS (in bp)", boxwex=0.9,
        at=c(1,2), outline=F, xlim=c(0.5,4.5), ylim=c(0,250), main="", cex.axis=1, font.axis=2)
boxplot(mrna_data_low$branchpt_len, mrna_data_high$branchpt_len, add= TRUE,col = c("#EB8A3E","#CB0000"), boxwex=0.9, at=c(3,4), outline=F, axes=F)
box(col="black", lwd=2)
dev.off()

lncrna_data_low$type <- "lncRNA (low TC)"
lncrna_data_high$type <- "lncRNA (high TC)"
mrna_data_low$type <- "mRNA (low TC)"
mrna_data_high$type <- "mRNA (high TC)"

library(tidyverse)
library(tidyquant)
library(ggdist)
library(ggthemes)

total_data <- data.frame(rbind(lncrna_data_low,lncrna_data_high,mrna_data_low,mrna_data_high))

total_data %>% filter(type %in% c("lncRNA (low TC)", "lncRNA (high TC)",
                                  "mRNA (low TC)", "mRNA (high TC)")) %>%
  ggplot(aes(x = factor(type), y = branchpt_len, fill = factor(type))) + 
  ggdist::stat_halfeye(adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA) +
  geom_boxplot(width=0.12,outlier.color = NA,alpha=0.5) +
  ggdist::stat_dots(side="left",justification = 1.1,binwidth=0.25) +
  scale_fill_tq() + theme_tq() +
  labs(title = "RainCloud Plot",x = "branchpoint length",y = "Highway Fuel",
    fill = "Cylinders") + coord_flip()