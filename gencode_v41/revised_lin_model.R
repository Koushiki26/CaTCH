setwd("/home/dell/Documents/koushiki/project/")

library(caret)
library(leaps)
library(car)
library(dplyr)
library(Metrics)
rsquare <- function(actual_y, predict_y){
  cor(actual_y, predict_y)^2
}

#All genes
main_data <- read.csv("data15.csv", header=T)
head(main_data)
main_data$TC <- main_data$Num_of_Transcripts/main_data$Num_of_Exons
main_data <- main_data[,-c(1:6,10,13,14,17:20)]
main_data <- na.omit(main_data)
names(main_data)[8] <- "chr_sig"
summary(main_data)
regsubset_out <- regsubsets(TC ~., data = main_data, nbest = 1, nvmax = NULL, 
                            force.in = NULL, force.out = NULL, method = "exhaustive")
regsubset_out
summary_out <- summary(regsubset_out)
as.data.frame(summary_out$outmat)
plot(regsubset_out, scale = "adjr2", main = "Adjusted R^2")

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

for(i in 1:1000){
  main_data$id <- 1:nrow(main_data)
  train_data <- main_data %>% dplyr::sample_frac(0.70)
  test_data <- dplyr::anti_join(main_data, train_data, by = 'id')
  test_data$TC1 <- ifelse(test_data$TC>0.3333333, 1, 0)
  lin_model1 <- lm(TC~Exon_Length, train_data)
  lin_predict1 <- predict(lin_model1,test_data[-9])
  predict1 <- ifelse(lin_predict1>0.3333333, 1, 0)
  auc_result[i,1] <- auc(test_data$TC1, predict1)
  lin_model2 <- lm(TC~Intron_Length, train_data)
  lin_predict2 <- predict(lin_model2,test_data[-9])
  predict2 <- ifelse(lin_predict2>0.3333333, 1, 0)
  auc_result[i,2] <- auc(test_data$TC1, predict2)
  lin_model3 <- lm(TC~codingpotential, train_data)
  lin_predict3 <- predict(lin_model3,test_data[-9])
  predict3 <- ifelse(lin_predict3>0.3333333, 1, 0)
  auc_result[i,3] <- auc(test_data$TC1, predict3)
  lin_model4 <- lm(TC~splice5, train_data)
  lin_predict4 <- predict(lin_model4,test_data[-9])
  predict4 <- ifelse(lin_predict4>0.3333333, 1, 0)
  auc_result[i,4] <- auc(test_data$TC1, predict4)
  lin_model5 <- lm(TC~chr_sig, train_data)
  lin_predict5 <- predict(lin_model5,test_data[-9])
  predict5 <- ifelse(lin_predict5>0.3333333, 1, 0)
  auc_result[i,5] <- auc(test_data$TC1, predict5)
  lin_model6 <- lm(TC~branchpt_len, train_data)
  lin_predict6 <- predict(lin_model6,test_data[-9])
  predict6 <- ifelse(lin_predict6>0.3333333, 1, 0)
  auc_result[i,6] <- auc(test_data$TC1, predict6)
  lin_model7 <- lm(TC~ssstrength5, train_data)
  lin_predict7 <- predict(lin_model7,test_data[-9])
  predict7<- ifelse(lin_predict7>0.3333333, 1, 0)
  auc_result[i,7] <- auc(test_data$TC1, predict7)
  lin_model8 <- lm(TC~ssstrength3, train_data)
  lin_predict8 <- predict(lin_model8,test_data[-9])
  predict8 <- ifelse(lin_predict8>0.3333333, 1, 0)
  auc_result[i,8] <- auc(test_data$TC1, predict8)
  lin_model9 <- lm(TC~Exon_Length+codingpotential, train_data)
  lin_predict9 <- predict(lin_model9,test_data[-9])
  predict9 <- ifelse(lin_predict9>0.3333333, 1, 0)
  auc_result[i,9] <- auc(test_data$TC1, predict9)
  lin_model10 <- lm(TC~Exon_Length+splice5, train_data)
  lin_predict10 <- predict(lin_model10,test_data[-9])
  predict10 <- ifelse(lin_predict10>0.3333333, 1, 0)
  auc_result[i,10] <- auc(test_data$TC1, predict10)
  lin_model11 <- lm(TC~Exon_Length+chr_sig, train_data)
  lin_predict11 <- predict(lin_model11,test_data[-9])
  predict11 <- ifelse(lin_predict11>0.3333333, 1, 0)
  auc_result[i,11] <- auc(test_data$TC1, predict11)
  lin_model12 <- lm(TC~Exon_Length+branchpt_len, train_data)
  lin_predict12 <- predict(lin_model12,test_data[-9])
  predict12 <- ifelse(lin_predict12>0.3333333, 1, 0)
  auc_result[i,12] <- auc(test_data$TC1, predict12)
  lin_model13 <- lm(TC~Exon_Length+Intron_Length, train_data)
  lin_predict13 <- predict(lin_model13,test_data[-9])
  predict13 <- ifelse(lin_predict13>0.3333333, 1, 0)
  auc_result[i,13] <- auc(test_data$TC1, predict13)
  lin_model14 <- lm(TC~codingpotential+chr_sig, train_data)
  lin_predict14 <- predict(lin_model14,test_data[-9])
  predict14 <- ifelse(lin_predict14>0.3333333, 1, 0)
  auc_result[i,14] <- auc(test_data$TC1, predict14)
  lin_model15 <- lm(TC~codingpotential+branchpt_len, train_data)
  lin_predict15 <- predict(lin_model15,test_data[-9])
  predict15 <- ifelse(lin_predict15>0.3333333, 1, 0)
  auc_result[i,15] <- auc(test_data$TC1, predict15)
  lin_model16 <- lm(TC~splice5+chr_sig, train_data)
  lin_predict16 <- predict(lin_model16,test_data[-9])
  predict16 <- ifelse(lin_predict16>0.3333333, 1, 0)
  auc_result[i,16] <- auc(test_data$TC1, predict16)
  lin_model17 <- lm(TC~codingpotential+splice5, train_data)
  lin_predict17 <- predict(lin_model17,test_data[-9])
  predict17 <- ifelse(lin_predict17>0.3333333, 1, 0)
  auc_result[i,17] <- auc(test_data$TC1, predict17)
  lin_model18 <- lm(TC~splice5+Intron_Length, train_data)
  lin_predict18 <- predict(lin_model18,test_data[-9])
  predict18 <- ifelse(lin_predict18>0.3333333, 1, 0)
  auc_result[i,18] <- auc(test_data$TC1, predict18)
  lin_model19 <- lm(TC~splice5+branchpt_len, train_data)
  lin_predict19 <- predict(lin_model19,test_data[-9])
  predict19 <- ifelse(lin_predict19>0.3333333, 1, 0)
  auc_result[i,19] <- auc(test_data$TC1, predict19)
  lin_model20 <- lm(TC~chr_sig+Intron_Length, train_data)
  lin_predict20 <- predict(lin_model20,test_data[-9])
  predict20 <- ifelse(lin_predict20>0.3333333, 1, 0)
  auc_result[i,20] <- auc(test_data$TC1, predict20)
  lin_model21 <- lm(TC~chr_sig+branchpt_len, train_data)
  lin_predict21 <- predict(lin_model21,test_data[-9])
  predict21 <- ifelse(lin_predict21>0.3333333, 1, 0)
  auc_result[i,21] <- auc(test_data$TC1, predict21)
  lin_model22 <- lm(TC~branchpt_len+Intron_Length, train_data)
  lin_predict22 <- predict(lin_model22,test_data[-9])
  predict22 <- ifelse(lin_predict22>0.3333333, 1, 0)
  auc_result[i,22] <- auc(test_data$TC1, predict22)
  lin_model23 <- lm(TC~Intron_Length+codingpotential, train_data)
  lin_predict23 <- predict(lin_model23,test_data[-9])
  predict23 <- ifelse(lin_predict23>0.3333333, 1, 0)
  auc_result[i,23] <- auc(test_data$TC1, predict23)
  lin_model24 <- lm(TC~Exon_Length+codingpotential+chr_sig, train_data)
  lin_predict24 <- predict(lin_model24,test_data[-9])
  predict24 <- ifelse(lin_predict24>0.3333333, 1, 0)
  auc_result[i,24] <- auc(test_data$TC1, predict24)
  lin_model25 <- lm(TC~Exon_Length+codingpotential+branchpt_len, train_data)
  lin_predict25 <- predict(lin_model25,test_data[-9])
  predict25 <- ifelse(lin_predict25>0.3333333, 1, 0)
  auc_result[i,25] <- auc(test_data$TC1, predict25)
  lin_model26 <- lm(TC~Exon_Length+splice5+chr_sig, train_data)
  lin_predict26 <- predict(lin_model26,test_data[-9])
  predict26 <- ifelse(lin_predict26>0.3333333, 1, 0)
  auc_result[i,26] <- auc(test_data$TC1, predict26)
  lin_model27 <- lm(TC~Exon_Length+codingpotential+splice5, train_data)
  lin_predict27 <- predict(lin_model27,test_data[-9])
  predict27 <- ifelse(lin_predict27>0.3333333, 1, 0)
  auc_result[i,27] <- auc(test_data$TC1, predict27)
  lin_model28 <- lm(TC~Exon_Length+splice5+Intron_Length, train_data)
  lin_predict28 <- predict(lin_model28,test_data[-9])
  predict28 <- ifelse(lin_predict28>0.3333333, 1, 0)
  auc_result[i,28] <- auc(test_data$TC1, predict28)
  lin_model29 <- lm(TC~Exon_Length+codingpotential+splice5+chr_sig, train_data)
  lin_predict29 <- predict(lin_model29,test_data[-9])
  predict29 <- ifelse(lin_predict29>0.3333333, 1, 0)
  auc_result[i,29] <- auc(test_data$TC1, predict29)
  lin_model30 <- lm(TC~Exon_Length+codingpotential+splice5+branchpt_len, train_data)
  lin_predict30 <- predict(lin_model30,test_data[-9])
  predict30 <- ifelse(lin_predict30>0.3333333, 1, 0)
  auc_result[i,30] <- auc(test_data$TC1, predict30)
  lin_model31 <- lm(TC~Exon_Length+splice5+chr_sig+Intron_Length, train_data)
  lin_predict31 <- predict(lin_model31,test_data[-9])
  predict31 <- ifelse(lin_predict31>0.3333333, 1, 0)
  auc_result[i,31] <- auc(test_data$TC1, predict31)
  lin_model32 <- lm(TC~Exon_Length+codingpotential+splice5+branchpt_len+chr_sig, train_data)
  lin_predict32 <- predict(lin_model32,test_data[-9])
  predict32 <- ifelse(lin_predict32>0.3333333, 1, 0)
  auc_result[i,32] <- auc(test_data$TC1, predict32)
  lin_model33 <- lm(TC~Intron_Length+Exon_Length+splice5+branchpt_len+chr_sig, train_data)
  lin_predict33 <- predict(lin_model33,test_data[-9])
  predict33 <- ifelse(lin_predict33>0.3333333, 1, 0)
  auc_result[i,33] <- auc(test_data$TC1, predict33)
  lin_model34 <- lm(TC~Exon_Length+codingpotential+splice5+branchpt_len+chr_sig+ssstrength5, train_data)
  lin_predict34 <- predict(lin_model34,test_data[-9])
  predict34 <- ifelse(lin_predict34>0.3333333, 1, 0)
  auc_result[i,34] <- auc(test_data$TC1, predict34)
  lin_model35 <- lm(TC~Exon_Length+codingpotential+splice5+branchpt_len+chr_sig+Intron_Length, train_data)
  lin_predict35 <- predict(lin_model35,test_data[-9])
  predict35 <- ifelse(lin_predict35>0.3333333, 1, 0)
  auc_result[i,35] <- auc(test_data$TC1, predict35)
  lin_model36 <- lm(TC~.-ssstrength3, train_data)
  lin_predict36 <- predict(lin_model36,test_data[-9])
  predict36 <- ifelse(lin_predict36>0.3333333, 1, 0)
  auc_result[i,36] <- auc(test_data$TC1, predict36)
  lin_model37 <- lm(TC~., train_data)
  lin_predict37 <- predict(lin_model37,test_data[-9])
  predict37 <- ifelse(lin_predict37>0.3333333, 1, 0)
  auc_result[i,37] <- auc(test_data$TC1, predict37)
}
summary(auc_result)
tiff("all_new1_lin.tiff",width = 3600, height = 1200, res = 300)
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
        outline=F, xlim=c(1,37), ylim=c(0.5,0.78), main="", cex.axis=1.6, font.axix=2)
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