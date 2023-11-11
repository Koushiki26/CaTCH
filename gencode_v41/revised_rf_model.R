#setwd("/home/dell/Documents/koushiki/project/")

library(randomForest)
library(dplyr)
library(pROC)
library(ggplot2)
library(Metrics)

main_data <- read.csv("data15.csv", header=T)
head(main_data)
main_data$TC <- main_data$Num_of_Transcripts/main_data$Num_of_Exons
main_data <- main_data[,-c(1:6,10,13,14,17:20)]
main_data <- na.omit(main_data)
names(main_data)[8] <- "chr_sig"
summary(main_data)

columns <- c("model1","model2","model3","model4","model5",
             "model6","model7","model8","model9","model10",
             "model11","model12","model13","model14","model15",
             "model16","model17","model18","model19","model20",
             "model21","model22","model23","model24","model25",
             "model26","model27","model28","model29","model30",
             "model31","model32","model33","model34","model35",
             "model36","model37")
auc_result <- data.frame(matrix(nrow=0, ncol = length(columns)))
colnames(auc_result) <- columns

for(i in 1:500){
  main_data$id <- 1:nrow(main_data)
  train <- main_data %>% dplyr::sample_frac(0.70)
  test <- dplyr::anti_join(main_data, train, by = 'id')
  test$TC1 <- ifelse(test$TC>0.3333333, 1, 0)
  model1 <- randomForest(TC~Exon_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict1TC <- predict(model1, newdata=test)
  predict1 <- ifelse(predict1TC>0.3333333, 1, 0)
  auc_result[i,1] <- auc(test$TC1, predict1)
  model2 <- randomForest(TC~codingpotential, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict2TC <- predict(model2, newdata=test)
  predict2 <- ifelse(predict2TC>0.3333333, 1, 0)
  auc_result[i,2] <- auc(test$TC1, predict2)
  model3 <- randomForest(TC~chr_sig, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict3TC <- predict(model3, newdata=test)
  predict3 <- ifelse(predict3TC>0.3333333, 1, 0)
  auc_result[i,3] <- auc(test$TC1, predict3)
  model4 <- randomForest(TC~Intron_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict4TC <- predict(model4, newdata=test)
  predict4 <- ifelse(predict4TC>0.3333333, 1, 0)
  auc_result[i,4] <- auc(test$TC1, predict4)
  model5 <- randomForest(TC~ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict5TC <- predict(model5, newdata=test)
  predict5 <- ifelse(predict5TC>0.3333333, 1, 0)
  auc_result[i,5] <- auc(test$TC1, predict5)
  model6 <- randomForest(TC~branchpt_len, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict6TC <- predict(model6, newdata=test)
  predict6 <- ifelse(predict6TC>0.3333333, 1, 0)
  auc_result[i,6] <- auc(test$TC1, predict6)
  model7 <- randomForest(TC~splice5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict7TC <- predict(model7, newdata=test)
  predict7 <- ifelse(predict7TC>0.3333333, 1, 0)
  auc_result[i,7] <- auc(test$TC1, predict7)
  model8 <- randomForest(TC~ssstrength3, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict8TC <- predict(model8, newdata=test)
  predict8 <- ifelse(predict8TC>0.3333333, 1, 0)
  auc_result[i,8] <- auc(test$TC1, predict8)
  model9 <- randomForest(TC~Exon_Length+codingpotential, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict9TC <- predict(model9, newdata=test)
  predict9 <- ifelse(predict9TC>0.3333333, 1, 0)
  auc_result[i,9] <- auc(test$TC1, predict9)
  model10 <- randomForest(TC~Exon_Length+chr_sig, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict10TC <- predict(model10, newdata=test)
  predict10 <- ifelse(predict10TC>0.3333333, 1, 0)
  auc_result[i,10] <- auc(test$TC1, predict10)
  model11 <- randomForest(TC~Exon_Length+Intron_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict11TC <- predict(model11, newdata=test)
  predict11 <- ifelse(predict11TC>0.3333333, 1, 0)
  auc_result[i,11] <- auc(test$TC1, predict11)
  model12 <- randomForest(TC~Exon_Length+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict12TC <- predict(model12, newdata=test)
  predict12 <- ifelse(predict12TC>0.3333333, 1, 0)
  auc_result[i,12] <- auc(test$TC1, predict12)
  model13 <- randomForest(TC~codingpotential+chr_sig, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict13TC <- predict(model13, newdata=test)
  predict13 <- ifelse(predict13TC>0.3333333, 1, 0)
  auc_result[i,13] <- auc(test$TC1, predict13)
  model14 <- randomForest(TC~codingpotential+Intron_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict14TC <- predict(model14, newdata=test)
  predict14 <- ifelse(predict14TC>0.3333333, 1, 0)
  auc_result[i,14] <- auc(test$TC1, predict14)
  model15 <- randomForest(TC~codingpotential+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict15TC <- predict(model15, newdata=test)
  predict15 <- ifelse(predict15TC>0.3333333, 1, 0)
  auc_result[i,15] <- auc(test$TC1, predict15)
  model16 <- randomForest(TC~chr_sig+Intron_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict16TC <- predict(model16, newdata=test)
  predict16 <- ifelse(predict16TC>0.3333333, 1, 0)
  auc_result[i,16] <- auc(test$TC1, predict16)
  model17 <- randomForest(TC~chr_sig+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict17TC <- predict(model17, newdata=test)
  predict17 <- ifelse(predict17TC>0.3333333, 1, 0)
  auc_result[i,17] <- auc(test$TC1, predict17)
  model18 <- randomForest(TC~Intron_Length+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict18TC <- predict(model18, newdata=test)
  predict18 <- ifelse(predict18TC>0.3333333, 1, 0)
  auc_result[i,18] <- auc(test$TC1, predict18)
  model19 <- randomForest(TC~Exon_Length+codingpotential+chr_sig, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict19TC <- predict(model19, newdata=test)
  predict19 <- ifelse(predict19TC>0.3333333, 1, 0)
  auc_result[i,19] <- auc(test$TC1, predict19)
  model20 <- randomForest(TC~Exon_Length+codingpotential+Intron_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict20TC <- predict(model20, newdata=test)
  predict20 <- ifelse(predict20TC>0.3333333, 1, 0)
  auc_result[i,20] <- auc(test$TC1, predict20)
  model21 <- randomForest(TC~Exon_Length+codingpotential+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict21TC <- predict(model21, newdata=test)
  predict21 <- ifelse(predict21TC>0.3333333, 1, 0)
  auc_result[i,21] <- auc(test$TC1, predict21)
  model22 <- randomForest(TC~Exon_Length+chr_sig+Intron_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict22TC <- predict(model22, newdata=test)
  predict22 <- ifelse(predict22TC>0.3333333, 1, 0)
  auc_result[i,22] <- auc(test$TC1, predict22)
  model23 <- randomForest(TC~Exon_Length+chr_sig+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict23TC <- predict(model23, newdata=test)
  predict23 <- ifelse(predict23TC>0.3333333, 1, 0)
  auc_result[i,23] <- auc(test$TC1, predict23)
  model24 <- randomForest(TC~Exon_Length+Intron_Length+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict24TC <- predict(model24, newdata=test)
  predict24 <- ifelse(predict24TC>0.3333333, 1, 0)
  auc_result[i,24] <- auc(test$TC1, predict24)
  model25 <- randomForest(TC~codingpotential+chr_sig+Intron_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict25TC <- predict(model25, newdata=test)
  predict25 <- ifelse(predict25TC>0.3333333, 1, 0)
  auc_result[i,25] <- auc(test$TC1, predict25)
  model26 <- randomForest(TC~codingpotential+chr_sig+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict26TC <- predict(model26, newdata=test)
  predict26 <- ifelse(predict26TC>0.3333333, 1, 0)
  auc_result[i,26] <- auc(test$TC1, predict26)
  model27 <- randomForest(TC~codingpotential+Intron_Length+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict27TC <- predict(model27, newdata=test)
  predict27 <- ifelse(predict27TC>0.3333333, 1, 0)
  auc_result[i,27] <- auc(test$TC1, predict27)
  model28 <- randomForest(TC~chr_sig+Intron_Length+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict28TC <- predict(model28, newdata=test)
  predict28 <- ifelse(predict28TC>0.3333333, 1, 0)
  auc_result[i,28] <- auc(test$TC1, predict28)
  model29 <- randomForest(TC~Exon_Length+codingpotential+chr_sig+Intron_Length, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict29TC <- predict(model29, newdata=test)
  predict29 <- ifelse(predict29TC>0.3333333, 1, 0)
  auc_result[i,29] <- auc(test$TC1, predict29)
  model30 <- randomForest(TC~Exon_Length+codingpotential+chr_sig+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict30TC <- predict(model30, newdata=test)
  predict30 <- ifelse(predict30TC>0.3333333, 1, 0)
  auc_result[i,30] <- auc(test$TC1, predict30)
  model31 <- randomForest(TC~Exon_Length+codingpotential+Intron_Length+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict31TC <- predict(model31, newdata=test)
  predict31 <- ifelse(predict31TC>0.3333333, 1, 0)
  auc_result[i,31] <- auc(test$TC1, predict31)
  model32 <- randomForest(TC~Exon_Length+chr_sig+Intron_Length+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict32TC <- predict(model32, newdata=test)
  predict32 <- ifelse(predict32TC>0.3333333, 1, 0)
  auc_result[i,32] <- auc(test$TC1, predict32)
  model33 <- randomForest(TC~codingpotential+chr_sig+Intron_Length+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict33TC <- predict(model33, newdata=test)
  predict33 <- ifelse(predict33TC>0.3333333, 1, 0)
  auc_result[i,33] <- auc(test$TC1, predict33)
  model34 <- randomForest(TC~Exon_Length+codingpotential+chr_sig+Intron_Length+ssstrength5, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict34TC <- predict(model34, newdata=test)
  predict34 <- ifelse(predict34TC>0.3333333, 1, 0)
  auc_result[i,34] <- auc(test$TC1, predict34)
  model35 <- randomForest(TC~Exon_Length+codingpotential+chr_sig+Intron_Length+ssstrength5+branchpt_len, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict35TC <- predict(model35, newdata=test)
  predict35 <- ifelse(predict35TC>0.3333333, 1, 0)
  auc_result[i,35] <- auc(test$TC1, predict35)
  model36 <- randomForest(TC~.-ssstrength3, data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict36TC <- predict(model36, newdata=test)
  predict36 <- ifelse(predict36TC>0.3333333, 1, 0)
  auc_result[i,36] <- auc(test$TC1, predict36)
  model37 <- randomForest(TC~., data = train,  mtry=2, importance=T, ntrees = 500, na.action = na.omit, type="regression")
  predict37TC <- predict(model37, newdata=test)
  predict37 <- ifelse(predict37TC>0.3333333, 1, 0)
  auc_result[i,37] <- auc(test$TC1, predict37)
}
summary(auc_result)
tiff("all_new1_rf.tiff",width = 3600, height = 1200, res = 300)
par(mar=c(0.5,2.5,0.5,0.5),cex.axis=1, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(auc_result, col=c("#fcfdbf","#fceeb0","#fddea0","#fecf92","#febf84",
                          "#feb078","#fe9f6d","#fc9065","#fa7f5e","#f7705c",
                          "#f1605d","#e95462","#de4968","#d2426f","#c43c75",
                          "#b73779","#a8327d","#992d80","#8c2981","#7e2482",
                          "#721f81","#641a80","#57157e","#3b0f70","#443983",
                          "#3d4d8a","#355f8d","#2d708e","#27808e","#21918c",
                          "#1fa188","#2ab07f","#35b779","#44bf70","#54c568",
                          "#7ad151","#a5db36"), 
        xlab="",xaxt='n', ylab="", boxwex=0.9,
        outline=F, xlim=c(1,37), ylim=c(0,0.4), main="", cex.axis=1.6, font.axix=2)
dev.off()


#ROC curve
library(caret)
library(leaps)
library(car)
library(dplyr)
library(pROC)
rsquare <- function(actual_y, predict_y){
  cor(actual_y, predict_y)^2
}

main_data <- read.csv("data15.csv", header=T)
head(main_data)
main_data$TC <- main_data$Num_of_Transcripts/main_data$Num_of_Exons
main_data <- main_data[,-c(1:6,10,13,14,17:20)]
main_data <- na.omit(main_data)
names(main_data)[8] <- "chr_sig"
summary(main_data)

main_data$id <- 1:nrow(main_data)
train <- main_data %>% dplyr::sample_frac(0.70)
lin_model <- lm(TC~Exon_Length+codingpotential+splice5+chr_sig, train)
lin_model$coefficients

#30% of gencode data
test1  <- dplyr::anti_join(main_data, train, by = 'id')
test1$TC1 <- ifelse(test1$TC>0.3333333, 1, 0)
predict_data1 <- predict(lin_model, newdata=test1)
test1$pred_TC <- predict_data1
predict1 <- ifelse(predict_data1>0.3333333, 1, 0)
rmse(test1$TC, predict(lin_model, newdata=test1))
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
predict_data2 <- predict(lin_model, newdata=test2)
test2$pred_TC <- predict_data2
predict2 <- ifelse(predict_data2>0.3333333, 1, 0)
rmse(test2$TC, predict(lin_model, newdata=test2))
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
predict_data3 <- predict(lin_model, newdata=test3)
test3$pred_TC <- predict_data3
predict3 <- ifelse(predict_data3>0.3333333, 1, 0)
rmse(test3$TC, predict(lin_model, newdata=test3))
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
predict_data4 <- predict(lin_model, newdata=test4)
test4$pred_TC <- predict_data4
predict4 <- ifelse(predict_data4>0.3333333, 1, 0)
rmse(test4$TC, predict(lin_model, newdata=test4))
auc(test4$TC1, predict4)
test_roc4 <- roc(test4$TC1~predict4, plot=T, print.auc=T)


tiff("roc_lin_new.tiff", width=1200, height=1200, res=300)
par(mar=c(2.5,2.5,0.5,0.5),font.lab=2,cex.lab=1.3,font.axis=2,cex.axis=1)
roc(test1$TC1~predict1, plot=T, print.auc=F, col="orangered2",lwd=2.5,xlim=c(1,0),ylim=c(0,1))
roc(test2$TC1~predict2, plot=T, print.auc=F, add=T, col="springgreen3",lwd=2.5)
roc(test3$TC1~predict3, plot=T, print.auc=F, add=T, col="blue4",lwd=2.5)
roc(test4$TC1~predict4, plot=T, print.auc=F, add=T, col="goldenrod1",lwd=2.5)
dev.off()


library(randomForest)
library(dplyr)
library(pROC)
library(ggplot2)
library(Metrics)

main_data <- read.csv("data15.csv", header=T)
head(main_data)
main_data$TC <- main_data$Num_of_Transcripts/main_data$Num_of_Exons
main_data <- main_data[,-c(1:6,10,13,14,17:20)]
main_data <- na.omit(main_data)
summary(main_data)
colnames(main_data)[8] <- "chr_sig"
data_rf <- randomForest(TC~., data = main_data, mtry=3, importance=T, na.action = na.omit, type="regression")
print(data_rf)

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

tiff("roc_rf_new.tiff", width=1200, height=1200, res=300)
par(mar=c(2.5,2.5,0.5,0.5),font.lab=2,cex.lab=1.3,font.axis=2,cex.axis=1)
roc(test1$TC1~predict1, plot=T, print.auc=F, col="orangered2",lwd=2.5,xlim=c(1,0),ylim=c(0,1))
roc(test2$TC1~predict2, plot=T, print.auc=F, add=T, col="springgreen3",lwd=2.5)
roc(test3$TC1~predict3, plot=T, print.auc=F, add=T, col="blue4",lwd=2.5)
roc(test4$TC1~predict4, plot=T, print.auc=F, add=T, col="goldenrod1",lwd=2.5)
dev.off()