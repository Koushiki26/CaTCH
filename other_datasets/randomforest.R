setwd("/home/dell/Documents/koushiki/project/")

library(randomForest)
library(dplyr)
library(pROC)
library(ggplot2)
library(Metrics)

main_data <- read.csv("data14.csv", header=T)
head(main_data)
main_data$TC <- main_data$Num_of_Transcripts/main_data$Num_of_Exons
main_data <- main_data[,-c(1:6,10,13,14,17:19)]
main_data <- na.omit(main_data)
summary(main_data)
data_rf <- randomForest(TC~., data = main_data, mtry=3, importance=T, na.action = na.omit, type="regression")
print(data_rf)

tiff("/home/dell/Documents/koushiki/project_new/trees.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5))
plot(data_rf, xlim=c(0,500), ylim=c(0.011,0.022),main="",cex.axis=1, 
     font.axis=2,col="#2988BC",lwd=2.5,cex.lab=1.3, font.lab=2)
dev.off()
Impdata <- as.data.frame(importance(data_rf))
Impdata$Var.Names <- row.names(Impdata)

tiff("/home/dell/Documents/koushiki/project_new/varimp.tiff",width = 1200, height = 1200, res = 300)
ggplot(Impdata, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="#7AA802") +
  geom_point(aes(size = IncNodePurity), color="#E1315B", alpha=0.6) +
  theme_light() + coord_flip() +xlab("") + ylab("%IncMSE") +
  theme(legend.position=c(0.8,0.8),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),legend.direction = "vertical",
    legend.box = "horizontal",
    axis.ticks.y = element_blank(),legend.text = element_text(size = 8, colour = "black", face="bold"),
    legend.title = element_text(size = 10, colour = "black", face="bold")) +
  theme(axis.title = element_text(size = 15, color = "black", face = "bold"),
        axis.text = element_text(size = 10, color = "black", face="bold"))
dev.off()
varImpPlot(data_rf)

main_data$id <- 1:nrow(main_data)
train <- main_data %>% dplyr::sample_frac(0.70)
model <- randomForest(TC~Exon_Length+codingpotential+chr_sig+Intron_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
print(model)

#30% of gencode data
test1  <- dplyr::anti_join(main_data, train, by = 'id')
test1$TC1 <- ifelse(test1$TC>0.3333333, 1, 0)
predict_data1 <- predict(model, newdata=test1)
test1$pred_TC <- predict_data1
predict1 <- ifelse(predict_data1>0.3333333, 1, 0)
rmse(test1$TC, predict(model, newdata=test1))
auc(test1$TC1, predict1)
test_roc1 <- roc(test1$TC1~predict1, plot=T, print.auc=T)

#gen_v38
test2 <- read.csv("../project_new/gencode/data9_gen.csv", header=T)
head(test2)
test2$TC <- test2$Num_of_Transcripts/test2$Num_of_Exons
test2 <- test2[,-c(1:5,12:15)]
test2 <- na.omit(test2)
summary(test2)
test2$TC1 <- ifelse(test2$TC>0.3333333, 1, 0)
predict_data2 <- predict(model, newdata=test2)
test2$pred_TC <- predict_data2
predict2 <- ifelse(predict_data2>0.3333333, 1, 0)
rmse(test2$TC, predict(model, newdata=test2))
auc(test2$TC1, predict2)
test_roc2 <- roc(test2$TC1~predict2, plot=T, print.auc=T)

#ensembl
test3 <- read.csv("../project_new/ensembl/data9_ens.csv", header=T)
head(test3)
test3$TC <- test3$Num_of_Transcripts/test3$Num_of_Exons
test3 <- test3[,-c(1:5,12:15)]
test3 <- na.omit(test3)
summary(test3)
test3$TC1 <- ifelse(test3$TC>0.3333333, 1, 0)
predict_data3 <- predict(model, newdata=test3)
test3$pred_TC <- predict_data3
predict3 <- ifelse(predict_data3>0.3333333, 1, 0)
rmse(test3$TC, predict(model, newdata=test3))
auc(test3$TC1, predict3)
test_roc3 <- roc(test3$TC1~predict3, plot=T, print.auc=T)

#t2t
test4 <- read.csv("../project_new/refseq_T2T/data8.csv", header=T)
head(test4)
test4$TC <- test4$Num_of_Transcripts/test4$Num_of_Exons
test4 <- test4[,-c(1:5,12:14)]
test4 <- na.omit(test4)
summary(test4)
test4$TC1 <- ifelse(test4$TC>0.3333333, 1, 0)
predict_data4 <- predict(model, newdata=test4)
test4$pred_TC <- predict_data4
predict4 <- ifelse(predict_data4>0.3333333, 1, 0)
rmse(test4$TC, predict(model, newdata=test4))
auc(test4$TC1, predict4)
test_roc4 <- roc(test4$TC1~predict4, plot=T, print.auc=T)

tiff("roc_rf.tiff", width=1200, height=1200, res=300)
par(mar=c(2.5,2.5,0.5,0.5),font.lab=2,cex.lab=1.3,font.axis=2,cex.axis=1)
roc(test1$TC1~predict1, plot=T, print.auc=F, col="orangered2",lwd=2.5,xlim=c(1,0),ylim=c(0,1))
roc(test2$TC1~predict2, plot=T, print.auc=F, add=T, col="springgreen3",lwd=2.5)
roc(test3$TC1~predict3, plot=T, print.auc=F, add=T, col="blue4",lwd=2.5)
roc(test4$TC1~predict4, plot=T, print.auc=F, add=T, col="goldenrod1",lwd=2.5)
dev.off()