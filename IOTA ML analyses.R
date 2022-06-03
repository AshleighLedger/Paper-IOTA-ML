#### IOTA ML paper analyses

#packages
library(nnet)
library(lme4)
library(glmnet)
library(brglm2)
library(caret)
library(auRoc)
library(car)
library(data.table)
library(boot)
library(metafor)
library(forestplot)
library(tiff)
library(xgboost)
library(ranger)
library(kernlab)
library(plyr)
library(dplyr)
library(plotrix)
library(readxl)
library(SDMTools)
library(bayesm)
library(VGAM)
library(ggplot2)
library(ggpubr)
library(ggformula)
library(doBy)


#### Without CA125 ####


mydata <- read.delim("iotadata_yang.txt")

factors <- c("outcome5", "outcome1", "loc10", "Shadows", "Ascites", "oncocenter", "colscore", "papnr", "Center")
for (i in factors) {
  mydata[, i] <- factor(mydata[, i])
}

#transformation
mydata$loglesdmax <- log(mydata$lesdmax)
mydata$logCA125 <- log(mydata$CA125)
mydata$papnrli <- as.numeric(mydata$papnr)
mydata$propsol2 <- I(mydata$propsol^2)

#normalize continuous variables
mydata$Agesd <- scale(mydata$Age, center=TRUE, scale=TRUE)
mydata$propsolsd <- scale(mydata$propsol, center=TRUE, scale=TRUE)
mydata$propsolsd2 <- I(mydata$propsolsd^2)
mydata$papnrlisd <- scale(mydata$papnrli, center=TRUE, scale=TRUE)
mydata$lesdmaxsd <- scale(mydata$lesdmax, center=TRUE, scale=TRUE)



#training data
train <- mydata[1:5914,]
train<-train[!(train$ADNEXexcluded==1),]

#validation dataset
test <- mydata[5915:10819,]
test <- test[!(test$operated120days==0),]


########################
#descriptive statistics#
########################

#stratified by outcome#
#nr of patients
by(traindf, traindf$outcome5, summary)
by(testdf, testdf$outcome5, summary)
#missing serum CA125
sum(is.na(train[train$outcome5==1, "CA125"]))
sum(is.na(test[test$outcome5==1, "CA125"]))
#stratified by dataset
summary(train)
summary(test)
#serum CA125
summary(traindf[, "ca125imp"])
summary(testdf[, "ca125imp"])

#####################
#Building the models#
#####################

#1. MLR

set.seed(123)
mm1 <- multinom(outcome5~oncocenter+Agesd+I(propsolsd^2)+I(propsolsd)+papnrlisd+loc10+Shadows+Ascites+loglesdmax, data=train)

#2. LINEAR MLR

mm2 <- multinom(outcome5~oncocenter+Agesd+propsolsd+papnr+loc10+Shadows+Ascites+lesdmaxsd, data=train)

#4. RIDGE MLR
#prepare data for glmnet
y_train_std <- data.matrix(train[,"outcome5"])
x_train_std <- data.matrix(train[, c("loglesdmax", "Agesd","oncocenter" ,"propsolsd","papnrlisd", "loc10", "Shadows", "Ascites", "propsolsd2")])
y_test_std <- data.matrix(test[,"outcome5"])
x_test_std <- data.matrix(test[, c("loglesdmax", "Agesd","oncocenter" ,"propsolsd","papnrlisd", "loc10", "Shadows", "Ascites", "propsolsd2")])

#cross-validation to tune lambda
set.seed(234)
cv_train_std <- cv.glmnet(x_train_std, y_train_std, nfolds=10, type.measure="auc", alpha=0, family="multinomial")
lambda <- cv_train_std$lambda.min
lambda
plot(cv_train_std)

#Build the model with optimal lambda value
mm4 <- glmnet(x_train_std, y_train_std, alpha=0, family="multinomial", lambda=lambda)

# 5. FIRTH MLR
mm5 <- pmlr(outcome5~oncocenter+Agesd+I(propsolsd)+I(propsolsd^2)+papnrlisd+loc10+Shadows+Ascites+loglesdmax, data=train,alpha = 0.05, penalized = TRUE, method = "likelihood")

# 7. RANDOM FOREST 
control <- trainControl(method="cv", number=10, classProbs=TRUE, summaryFunction=multiClassSummary, search="random")
set.seed(1234)
mm7 <- train(make.names(outcome5)~oncocenter+Age+propsol+papnr+loc10+Shadows+Ascites+lesdmax, 
             data = train, 
             method = "ranger",
             trControl = control,
             tuneLength=30,
             metric="logLoss")
mm7

# 8. XGBoost 
#PREPARE DATA FOR XGBOOST
train22 <- subset(train, select=c("outcome5", "oncocenter", "lesdmax", "papnr", "loc10", "Shadows", "Ascites", "Age", "propsol"))
test22 <-subset(test, select=c("outcome5", "oncocenter", "lesdmax", "papnr", "loc10", "Shadows", "Ascites", "Age", "propsol"))
dummies <- dummyVars(~., data=train22[,-1])
c2 <- predict(dummies, train22[, -1])
d_training <- as.data.frame(cbind(train22$outcome5, c2))
dummies <- dummyVars(~., data=test22[, -1])
c2 <- predict(dummies, test22[,-1])
d_test <- as.data.frame(cbind(test22$outcome5, c2))

x_train <- xgb.DMatrix(as.matrix(d_training %>% select (-V1)))
y_train <- as.factor(make.names(d_training$V1))

x_test <- xgb.DMatrix(as.matrix(d_test %>% select(-V1)))
y_test <- as.factor(make.names(d_test$V1))

xgb_trcontrol = trainControl(method = "cv", number = 10, allowParallel = TRUE, 
                             verboseIter = FALSE, returnData = FALSE, summaryFunction=multiClassSummary,   
                             classProbs=TRUE, search="random")

set.seed(045)
mm8 = train(x_train, y_train, trControl = xgb_trcontrol, 
            method = "xgbTree", tuneLength=30, metric="logLoss")
mm8


######################
#SETUP FOR NN AND SVM#
######################
train2 <- subset(train, select=c("outcome5", "oncocenter", "lesdmaxsd", "papnr", "loc10", "Shadows", "Ascites", "Agesd", "propsolsd"))
test2 <- subset(test, select=c("outcome5", "oncocenter", "lesdmaxsd", "papnr", "loc10", "Shadows", "Ascites", "Agesd", "propsolsd"))
#one hot encoding
dummies <- dummyVars(~., data=train2[,-1])
c2 <- predict(dummies, train2[, -1])
d_training <- as.data.frame(cbind(train2$outcome5, c2))
d_training$V1 <- as.factor(d_training$V1)
dummies <- dummyVars(~., data=test2[, -1])
c2 <- predict(dummies, test2[,-1])
d_test <- as.data.frame(cbind(test2$outcome5, c2))
d_test$V1 <- as.factor(d_test$V1)


# 9. NEURAL NETWORKS 
ctrl <- trainControl(method="cv",   # 10fold cross validation
                     number=10,         
                     summaryFunction=multiClassSummary,
                     search="random",
                     classProbs=TRUE)
set.seed(789)
mm9 <- 
  train(make.names(V1)~. ,
        data = d_training,
        method = "nnet",
        trControl = ctrl,
        trace = FALSE,
        maxit = 1000,
        linout = FALSE,
        tuneLength = 30,
        metric="logLoss") 
mm9


# 10. SVM 
set.seed(123)
system.time(mm10 <- train(form=make.names(V1)~.,
                          method = "svmRadial",# Radial basis function kernel
                          tuneLength = 30,
                          metric = "logLoss",
                          data=d_training,
                          trControl=ctrl) ) 
mm10



########################
##APPARENT PERFORMANCE##
########################

#Merge the smallest centers into 1
levels(train$Center)[levels(train$Center)=="FIT"] <-"Other"
levels(train$Center)[levels(train$Center)=="GIT"] <-"Other"
levels(train$Center)[levels(train$Center)=="OCA"] <-"Other"
levels(train$Center)[levels(train$Center)=="VIT"] <-"Other"


#AUC BENIGN VS MALIGNANT

# MLR
# Risk estimates
mm1_predapp <- predict(mm1,  type="probs")
train$mm1_predapp <- rowSums (mm1_predapp[,2:5])

apparentm1auc <- AUC.IOTA(pred = mm1_predapp, outcome=outcome1, center=Center, data=train)

#FOREST PLOT
apparentm1auc$Plot
apparentm1auc$Plot[24, "RRauc"] <- "   "
apparentm1auc$Plot[26, "RRauc"] <- "        (0.89 to 0.96)"
apparentm1auc$Plot[26, "RRprev"] <- "   "
apparentm1auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven", "Rome", "Malmö", "Genk", "Monza", "Prague", 
                                      "Bologna","Milan 1 ", "Lublin", "Cagliari","Milan 3", "Stockholm", 
                                      "London", "Naples","Paris","Lund", "Beijing","Other", "Maurepas","Udin","Barcelona", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")

tiff("mlr auc apparent.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(apparentm1auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentm1auc$dataPlot$AUC,
           lower = apparentm1auc$dataPlot$LL,
           upper = apparentm1auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentm1auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "multinomial LR: apparent performance",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentm1auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentm1auc$Performance$AUC)

dev.off()

#LR WITHOUT TRANSFORMATIONS
mm2_predapp <- predict(mm22,  type="probs")
train$mm2_predapp <- rowSums (mm2_predapp[,2:5])

apparentm2auc <- AUC.IOTA(pred = mm2_predapp, outcome=outcome1, center=Center, data=train)

#FOREST PLOT
apparentm2auc$Plot
apparentm2auc$Plot[24, "RRauc"] <- "   "
apparentm2auc$Plot[26, "RRauc"] <- "        (0.88 to 0.96)"
apparentm2auc$Plot[26, "RRprev"] <- "   "
apparentm2auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                      "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France",
                                      "Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("linear mlr auc apparent.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(apparentm2auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentm2auc$dataPlot$AUC,
           lower = apparentm2auc$dataPlot$LL,
           upper = apparentm2auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentm2auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "multinomial LR w/o transformations: apparent performance",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentm2auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentm2auc$Performance$AUC)

dev.off()


#ridge
mm4_predapp <- predict(mm4, x_train_std, type="response", s=lambda)
mm4_predapp <- data.frame(mm4_predapp)
#calculating overall risk of malignancy
train$mm4_predapp <- rowSums (mm4_predapp[,2:5])

apparentm4auc <- AUC.IOTA(pred = mm4_predapp, outcome=outcome1, center=Center, data=train)

#FOREST PLOT
apparentm4auc$Plot
apparentm4auc$Plot[24, "RRauc"] <- "   "
apparentm4auc$Plot[26, "RRauc"] <- "        (0.88 to 0.96)"
apparentm4auc$Plot[26, "RRprev"] <- "   "
apparentm4auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                      "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France",
                                      "Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("ridge auc apparent.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(apparentm4auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentm4auc$dataPlot$AUC,
           lower = apparentm4auc$dataPlot$LL,
           upper = apparentm4auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentm4auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "multinomial ridge LR: apparent performance",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentm4auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentm4auc$Performance$AUC)
dev.off()

#firth

mm5_predapp
train$mm5_predapp <- rowSums (mm5_predapp[,2:5])

apparentm5auc <- AUC.IOTA(pred = mm5_predapp, outcome=outcome1, center=Center, data=train)

#FOREST PLOT
apparentm5auc$Plot
apparentm5auc$Plot[24, "RRauc"] <- "   "
apparentm5auc$Plot[26, "RRauc"] <- "        (0.89 to 0.96)"
apparentm5auc$Plot[26, "RRprev"] <- "   "
apparentm5auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                      "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France",
                                      "Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("firth auc apparent.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(apparentm5auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentm5auc$dataPlot$AUC,
           lower = apparentm5auc$dataPlot$LL,
           upper = apparentm5auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentm5auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "multinomial firth LR: apparent performance",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentm5auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentm5auc$Performance$AUC)

dev.off()


#RF

mm7_predapp <- predict(mm7,  type="prob")
newpred <- predict(mm7, newdata=train, type="prob")
train$mm7_predapp <- rowSums (mm7_predapp[,2:5])

apparentm7auc <- AUC.IOTA(pred = mm7_predapp, outcome=outcome1, center=Center, data=train)

#FOREST PLOT
apparentm7auc$Plot
apparentm7auc$Plot[24, "RRauc"] <- "   "
apparentm7auc$Plot[26, "RRauc"] <- "        (0.95 to 0.98)"
apparentm7auc$Plot[26, "RRprev"] <- "   "
apparentm7auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                      "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France",
                                      "Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("rf auc apparent.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(apparentm7auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentm7auc$dataPlot$AUC,
           lower = apparentm7auc$dataPlot$LL,
           upper = apparentm7auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentm7auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "multinomial RF: apparent performance",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentm7auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentm7auc$Performance$AUC)

dev.off()

#XG
mm8_predapp <- predict(mm8, newdata = x_train, type="prob")
train$mm8_predapp <- rowSums (mm8_predapp[,2:5])


apparentm8auc <- AUC.IOTA(pred = mm8_predapp, outcome=train$outcome1, center=train$Center, data=train)

#FOREST PLOT
apparentm8auc$Plot
apparentm8auc$Plot[24, "RRauc"] <- "   "
apparentm8auc$Plot[26, "RRauc"] <- "        (0.90 to 0.97)"
apparentm8auc$Plot[26, "RRprev"] <- "   "
apparentm8auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                      "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France",
                                      "Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("xg auc apparent.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(apparentm8auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentm8auc$dataPlot$AUC,
           lower = apparentm8auc$dataPlot$LL,
           upper = apparentm8auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentm8auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "multinomial XGBoost: apparent performance",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.50,0.75, 1),
           clip = c(0.50, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentm8auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentm8auc$Performance$AUC)

dev.off()

#NN
mm9_predapp <-predict(object = mm9,newdata = d_training,type = "prob")
#calculation of overall risk of malignancy
train$mm9_predapp <- rowSums (mm9_predapp[,2:5])

apparentm9auc <- AUC.IOTA(pred = mm9_predapp, outcome=train$outcome1, center=train$Center, data=train)

#FOREST PLOT
apparentm9auc$Plot
apparentm9auc$Plot[24, "RRauc"] <- "   "
apparentm9auc$Plot[26, "RRauc"] <- "        (0.88 to 0.97)"
apparentm9auc$Plot[26, "RRprev"] <- "   "
apparentm9auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                      "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France",
                                      "Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("multi nn auc apparent.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(apparentm9auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentm9auc$dataPlot$AUC,
           lower = apparentm9auc$dataPlot$LL,
           upper = apparentm9auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentm9auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "multinomial NN: apparent performance",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentm9auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentm9auc$Performance$AUC)

dev.off()

#SVM
mm10_predapp <-predict(object = mm10,newdata = d_training %>% select(-V1), type = "prob")
train$mm10_predapp <- rowSums (mm10_predapp[,2:5])

apparentm10auc <- AUC.IOTA(pred = mm10_predapp, outcome=train$outcome1, center=train$Center, data=train)

#FOREST PLOT
apparentm10auc$Plot
apparentm10auc$Plot[24, "RRauc"] <- "   "
apparentm10auc$Plot[26, "RRauc"] <- "        (0.86 to 0.96)"
apparentm10auc$Plot[26, "RRprev"] <- "   "
apparentm10auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                       "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France",
                                       "Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("svm auc apparent 2.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(apparentm10auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentm10auc$dataPlot$AUC,
           lower = apparentm10auc$dataPlot$LL,
           upper = apparentm10auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentm10auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "multinomial SVM: apparent performance",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.50, 0.75,1),
           clip = c(0.50, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentm10auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentm10auc$Performance$AUC)

dev.off()


#Forest plot for each model
NA.forest <- apparentm1auc$Performance[1,]
NA.forest <- NA
Summary.AUC <- rbind(NA.forest, apparentm1auc$Performance[1,], apparentm4auc$Performance[1,], apparentm5auc$Performance[1,],apparentm2auc$Performance[1,], apparentm7auc$Performance[1,], apparentm8auc$Performance[1,], apparentm9auc$Performance[1,], apparentm10auc$Performance[1,])
Summary.AUC$Model <- c('', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'Random forest', "XGBoost", "Neural network", "Support vector machine")
Summary.AUC.PI <- rbind(NA.forest, apparentm1auc$Performance[2,], apparentm4auc$Performance[2,], apparentm5auc$Performance[2,], apparentm2auc$Performance[2,] , apparentm7auc$Performance[2,], apparentm8auc$Performance[2,], apparentm9auc$Performance[2,], apparentm10auc$Performance[2,])
Summary.AUC.PI$Model <- c('', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'Random forest', "XGBoost", "Neural network", "Support vector machine")
tabletext <- cbind(
  c('Model', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', "Random forest", "XGBoost",'Neural network', "Support vector machine"),
  c('AUROC (95% CI)', 
    paste(format(round(apparentm1auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentm1auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentm1auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentm4auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentm4auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentm4auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentm5auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentm5auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentm5auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentm2auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentm2auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentm2auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentm7auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentm7auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentm7auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentm8auc$Performance$AUC[1], 2), nsmall=2), " (", format(round(apparentm8auc$Performance$LL[1], 2), nsmall=2), " to ", format(round(apparentm8auc$Performance$UL[1], 2), nsmall=2), ")", sep= ""),
    paste(format(round(apparentm9auc$Performance$AUC[1], 2), nsmall=2), " (", format(round(apparentm9auc$Performance$LL[1], 2), nsmall=2), " to ", format(round(apparentm9auc$Performance$UL[1], 2), nsmall=2), ")", sep=""),
    paste(format(round(apparentm10auc$Performance$AUC[1], 2), nsmall=2), " (", format(round(apparentm10auc$Performance$LL[1], 2), nsmall=2), " to ", format(round(apparentm10auc$Performance$UL[1], 2), nsmall=2), ")", sep="")),
  c('95% PI', 
    paste0("(", format(round(apparentm1auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentm1auc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(apparentm4auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentm4auc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(apparentm5auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentm5auc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(apparentm2auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentm2auc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(apparentm7auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentm7auc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(apparentm8auc$Performance$LL[2], 2), nsmall=2), " to " , format(round(apparentm8auc$Performance$UL[2], 2), nsmall=2), ")"),
    paste0("(", format(round(apparentm9auc$Performance$LL[2], 2), nsmall=2), " to " , format(round(apparentm9auc$Performance$UL[2], 2), nsmall=2), ")"),
    paste0("(", format(round(apparentm10auc$Performance$LL[2], 2), nsmall=2), " to " , format(round(apparentm10auc$Performance$UL[2], 2), nsmall=2), ")")))


tiff("discrimination overall apparent 2.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE),
           xlab = "AUROC",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.80, 0.90, 1), xlog = TRUE, clip = c(0.80, 1))       
dev.off()



#PDI
#MLR
mm1pdiapp <- ests(y=train$outcome5, d=mm1_predapp, acc="pdi", level=0.95, method="prob", k=5)
mm1pdiapp

#LINEAR MLR
mm2pdiapp <- ests(y=train$outcome5, d=mm2_predapp, acc="pdi", level=0.95, method="prob", k=5)
mm2pdiapp


#RIDGE LOGISTIC REGRESSION
mm4pdiapp <- ests(y=train$outcome5, d=mm4_predapp, acc="pdi", level=0.95, method="prob", k=5)
mm4pdiapp

#FIRTH LOGISTIC REGRESSION
mm5pdiapp <- ests(y=train$outcome5, d=mm5_predapp, acc="pdi", level=0.95, method="prob", k=5)
mm5pdiapp


#RANDOM FOREST
mm7pdiapp <- ests(y=train$outcome5, d=mm7_predapp, acc="pdi", level=0.95, method="prob", k=5)
mm7pdiapp

#XGBOOST
mm8pdiapp <- ests(y=train$outcome5, d=mm8_predapp, acc="pdi", level=0.95, method="prob", k=5)
mm8pdiapp

#NEURAL NETWORK
mm9pdiapp <- ests(y=train$outcome5, d=mm9_predapp, acc="pdi", level=0.95, method="prob", k=5)
mm9pdiapp
#SVM
mm10pdiapp <- ests(y=train$outcome5, d=mm10_predapp, acc="pdi", level=0.95, method="prob", k=5)
mm10pdiapp

#pairwise AUC
#MLR
#Pair 1: benign vs borderline
train$bvsb <- mm1_predapp[, 1] / (mm1_predapp[,1] + mm1_predapp[,2])
Df = data.frame(p = train$bvsb, y = train$outcome5, stringsAsFactors = F)
AUCbvsb <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsb) <- c('AUC', 'LL', 'UL')
AUCbvsb <- data.frame(AUCbvsb)
AUCbvsblr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 2], method = "pepe")
AUCbvsb[1, 1:3] <- AUCbvsblr
AUCbvsb$logit.AUC <- logit(AUCbvsb$AUC)
AUCbvsb$logit.se  <- (logit(AUCbvsb$AUC) - logit(AUCbvsb$LL))/1.96
AUCbvsb$logit.var <- AUCbvsb$logit.se^2
AUCbvsb
#Pair 2: benign vs stage I
train$bvss1 <- mm1_predapp[,1] / (mm1_predapp[,1] +mm1_predapp[,3])
Df = data.frame(p = train$bvss1, y = train$outcome5, stringsAsFactors = F)
AUCbvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1) <- c('AUC', 'LL', 'UL')
AUCbvss1 <- data.frame(AUCbvss1)
AUCbvss1lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 3], method = "pepe")
AUCbvss1[1, 1:3] <- AUCbvss1lr
#Pair 3: benign vs stage II-IV
train$bvss2 <- mm1_predapp[,1] / (mm1_predapp[,1] +mm1_predapp[,4])
Df = data.frame(p = train$bvss2, y = train$outcome5, stringsAsFactors = F)
AUCbvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2) <- c('AUC', 'LL', 'UL')
AUCbvss2 <- data.frame(AUCbvss2)
AUCbvss2lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 4], method = "pepe")
AUCbvss2[1, 1:3] <- AUCbvss2lr
AUCbvss2
#Pair 4: benign vs secondary metastatic
train$bvssm <- mm1_predapp[,1] / (mm1_predapp[,1] +mm1_predapp[,5])
Df = data.frame(p = train$bvssm, y = train$outcome5, stringsAsFactors = F)
AUCbvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssm) <- c('AUC', 'LL', 'UL')
AUCbvssm <- data.frame(AUCbvssm)
AUCbvssmlr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 5], method = "pepe")
AUCbvssm[1, 1:3] <- AUCbvssmlr
AUCbvssm
#Pair 5: borderline vs stage I
train$bdvss1 <- mm1_predapp[,2] / (mm1_predapp[,2] +mm1_predapp[,3])
Df = data.frame(p = train$bdvss1, y = train$outcome5, stringsAsFactors = F)
AUCbdvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1) <- c('AUC', 'LL', 'UL')
AUCbdvss1 <- data.frame(AUCbdvss1)
AUCbdvss1lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 3], method = "pepe")
AUCbdvss1[1, 1:3] <- AUCbdvss1lr
AUCbdvss1
#Pair 6: borderline vs stage II-IV
train$bdvss2 <- mm1_predapp[,2] / (mm1_predapp[,2] +mm1_predapp[,4])
Df = data.frame(p = train$bdvss2, y = train$outcome5, stringsAsFactors = F)
AUCbdvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2) <- c('AUC', 'LL', 'UL')
AUCbdvss2 <- data.frame(AUCbdvss2)
AUCbdvss2lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 4], method = "pepe")
AUCbdvss2[1, 1:3] <- AUCbdvss2lr
AUCbdvss2
#pair 7: borderline vs secondary metastatic
train$bdvssm <- mm1_predapp[,2] / (mm1_predapp[,2] +mm1_predapp[,5])
Df = data.frame(p = train$bdvssm, y = train$outcome5, stringsAsFactors = F)
AUCbdvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssm) <- c('AUC', 'LL', 'UL')
AUCbdvssm <- data.frame(AUCbdvssm)
AUCbdvssmlr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 5], method = "pepe")
AUCbdvssm[1, 1:3] <- AUCbdvssmlr
AUCbdvssm
#pair 8: stage I vs stage II-IV
train$s1vss2 <- mm1_predapp[,3] / (mm1_predapp[,3] +mm1_predapp[,4])
Df = data.frame(p = train$s1vss2, y = train$outcome5, stringsAsFactors = F)
AUCs1vss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2) <- c('AUC', 'LL', 'UL')
AUCs1vss2 <- data.frame(AUCs1vss2)
AUCs1vss2lr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 4], method = "pepe")
AUCs1vss2[1, 1:3] <- AUCs1vss2lr
AUCs1vss2
#pair 9: stage I vs secondary metastatic
train$s1vssm <- mm1_predapp[,3] / (mm1_predapp[,3] +mm1_predapp[,5])
Df = data.frame(p = train$s1vssm, y = train$outcome5, stringsAsFactors = F)
AUCs1vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssm) <- c('AUC', 'LL', 'UL')
AUCs1vssm <- data.frame(AUCs1vssm)
AUCs1vssmlr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 5], method = "pepe")
AUCs1vssm[1, 1:3] <- AUCs1vssmlr
AUCs1vssm
#pair 10: stage II-IV vs secondary metastatic
train$s2vssm <- mm1_predapp[,4] / (mm1_predapp[,4] +mm1_predapp[,5])
Df = data.frame(p = train$s2vssm, y = train$outcome5, stringsAsFactors = F)
AUCs2vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssm) <- c('AUC', 'LL', 'UL')
AUCs2vssm <- data.frame(AUCs2vssm)
AUCs2vssmlr<- auc.nonpara.mw(Df$p[Df$y == 4], Df$p[Df$y == 5], method = "pepe")
AUCs2vssm[1, 1:3] <- AUCs2vssmlr
AUCs2vssm

#MLR
#Pair 1: benign vs borderline
train$bvsb <- mm2_predapp[, 1] / (mm2_predapp[,1] + mm2_predapp[,2])
Df = data.frame(p = train$bvsb, y = train$outcome5, stringsAsFactors = F)
AUCbvsb <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsb) <- c('AUC', 'LL', 'UL')
AUCbvsb <- data.frame(AUCbvsb)
AUCbvsblr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 2], method = "pepe")
AUCbvsb[1, 1:3] <- AUCbvsblr
#Pair 2: benign vs stage I
train$bvss1 <- mm2_predapp[,1] / (mm2_predapp[,1] +mm2_predapp[,3])
Df = data.frame(p = train$bvss1, y = train$outcome5, stringsAsFactors = F)
AUCbvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1) <- c('AUC', 'LL', 'UL')
AUCbvss1 <- data.frame(AUCbvss1)
AUCbvss1lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 3], method = "pepe")
AUCbvss1[1, 1:3] <- AUCbvss1lr
AUCbvss1$logit.AUC <- logit(AUCbvss1$AUC)
AUCbvss1$logit.se  <- (logit(AUCbvss1$AUC) - logit(AUCbvss1$LL))/1.96
AUCbvss1$logit.var <- AUCbvss1$logit.se^2
AUCbvss1
#Pair 3: benign vs stage II-IV
train$bvss2 <- mm2_predapp[,1] / (mm2_predapp[,1] +mm2_predapp[,4])
Df = data.frame(p = train$bvss2, y = train$outcome5, stringsAsFactors = F)
AUCbvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2) <- c('AUC', 'LL', 'UL')
AUCbvss2 <- data.frame(AUCbvss2)
AUCbvss2lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 4], method = "pepe")
AUCbvss2[1, 1:3] <- AUCbvss2lr
AUCbvss2
#Pair 4: benign vs secondary metastatic
train$bvssm <- mm2_predapp[,1] / (mm2_predapp[,1] +mm2_predapp[,5])
Df = data.frame(p = train$bvssm, y = train$outcome5, stringsAsFactors = F)
AUCbvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssm) <- c('AUC', 'LL', 'UL')
AUCbvssm <- data.frame(AUCbvssm)
AUCbvssmlr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 5], method = "pepe")
AUCbvssm[1, 1:3] <- AUCbvssmlr
AUCbvssm
#Pair 5: borderline vs stage I
train$bdvss1 <- mm2_predapp[,2] / (mm2_predapp[,2] +mm2_predapp[,3])
Df = data.frame(p = train$bdvss1, y = train$outcome5, stringsAsFactors = F)
AUCbdvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1) <- c('AUC', 'LL', 'UL')
AUCbdvss1 <- data.frame(AUCbdvss1)
AUCbdvss1lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 3], method = "pepe")
AUCbdvss1[1, 1:3] <- AUCbdvss1lr
AUCbdvss1
#Pair 6: borderline vs stage II-IV
train$bdvss2 <- mm2_predapp[,2] / (mm2_predapp[,2] +mm2_predapp[,4])
Df = data.frame(p = train$bdvss2, y = train$outcome5, stringsAsFactors = F)
AUCbdvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2) <- c('AUC', 'LL', 'UL')
AUCbdvss2 <- data.frame(AUCbdvss2)
AUCbdvss2lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 4], method = "pepe")
AUCbdvss2[1, 1:3] <- AUCbdvss2lr
AUCbdvss2
#pair 7: borderline vs secondary metastatic
train$bdvssm <- mm2_predapp[,2] / (mm2_predapp[,2] +mm2_predapp[,5])
Df = data.frame(p = train$bdvssm, y = train$outcome5, stringsAsFactors = F)
AUCbdvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssm) <- c('AUC', 'LL', 'UL')
AUCbdvssm <- data.frame(AUCbdvssm)
AUCbdvssmlr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 5], method = "pepe")
AUCbdvssm[1, 1:3] <- AUCbdvssmlr
AUCbdvssm
#pair 8: stage I vs stage II-IV
train$s1vss2 <- mm2_predapp[,3] / (mm2_predapp[,3] +mm2_predapp[,4])
Df = data.frame(p = train$s1vss2, y = train$outcome5, stringsAsFactors = F)
AUCs1vss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2) <- c('AUC', 'LL', 'UL')
AUCs1vss2 <- data.frame(AUCs1vss2)
AUCs1vss2lr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 4], method = "pepe")
AUCs1vss2[1, 1:3] <- AUCs1vss2lr
AUCs1vss2
#pair 9: stage I vs secondary metastatic
train$s1vssm <- mm2_predapp[,3] / (mm2_predapp[,3] +mm2_predapp[,5])
Df = data.frame(p = train$s1vssm, y = train$outcome5, stringsAsFactors = F)
AUCs1vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssm) <- c('AUC', 'LL', 'UL')
AUCs1vssm <- data.frame(AUCs1vssm)
AUCs1vssmlr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 5], method = "pepe")
AUCs1vssm[1, 1:3] <- AUCs1vssmlr
AUCs1vssm
#pair 10: stage II-IV vs secondary metastatic
train$s2vssm <- mm2_predapp[,4] / (mm2_predapp[,4] +mm2_predapp[,5])
Df = data.frame(p = train$s2vssm, y = train$outcome5, stringsAsFactors = F)
AUCs2vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssm) <- c('AUC', 'LL', 'UL')
AUCs2vssm <- data.frame(AUCs2vssm)
AUCs2vssmlr<- auc.nonpara.mw(Df$p[Df$y == 4], Df$p[Df$y == 5], method = "pepe")
AUCs2vssm[1, 1:3] <- AUCs2vssmlr
AUCs2vssm


#RIDGE
#Pair 1: benign vs borderline
train$bvsb <- mm4_predapp[, 1] / (mm4_predapp[,1] + mm4_predapp[,2])
Df = data.frame(p = train$bvsb, y = train$outcome5, stringsAsFactors = F)
AUCbvsb <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsb) <- c('AUC', 'LL', 'UL')
AUCbvsb <- data.frame(AUCbvsb)
AUCbvsblr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 2], method = "pepe")
AUCbvsb[1, 1:3] <- AUCbvsblr
#Pair 2: benign vs stage I
train$bvss1 <- mm4_predapp[,1] / (mm4_predapp[,1] +mm4_predapp[,3])
Df = data.frame(p = train$bvss1, y = train$outcome5, stringsAsFactors = F)
AUCbvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1) <- c('AUC', 'LL', 'UL')
AUCbvss1 <- data.frame(AUCbvss1)
AUCbvss1lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 3], method = "pepe")
AUCbvss1[1, 1:3] <- AUCbvss1lr
AUCbvss1$logit.AUC <- logit(AUCbvss1$AUC)
AUCbvss1$logit.se  <- (logit(AUCbvss1$AUC) - logit(AUCbvss1$LL))/1.96
AUCbvss1$logit.var <- AUCbvss1$logit.se^2
AUCbvss1
#Pair 3: benign vs stage II-IV
train$bvss2 <- mm4_predapp[,1] / (mm4_predapp[,1] +mm4_predapp[,4])
Df = data.frame(p = train$bvss2, y = train$outcome5, stringsAsFactors = F)
AUCbvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2) <- c('AUC', 'LL', 'UL')
AUCbvss2 <- data.frame(AUCbvss2)
AUCbvss2lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 4], method = "pepe")
AUCbvss2[1, 1:3] <- AUCbvss2lr
AUCbvss2
#Pair 4: benign vs secondary metastatic
train$bvssm <- mm4_predapp[,1] / (mm4_predapp[,1] +mm4_predapp[,5])
Df = data.frame(p = train$bvssm, y = train$outcome5, stringsAsFactors = F)
AUCbvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssm) <- c('AUC', 'LL', 'UL')
AUCbvssm <- data.frame(AUCbvssm)
AUCbvssmlr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 5], method = "pepe")
AUCbvssm[1, 1:3] <- AUCbvssmlr
AUCbvssm
#Pair 5: borderline vs stage I
train$bdvss1 <- mm4_predapp[,2] / (mm4_predapp[,2] +mm4_predapp[,3])
Df = data.frame(p = train$bdvss1, y = train$outcome5, stringsAsFactors = F)
AUCbdvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1) <- c('AUC', 'LL', 'UL')
AUCbdvss1 <- data.frame(AUCbdvss1)
AUCbdvss1lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 3], method = "pepe")
AUCbdvss1[1, 1:3] <- AUCbdvss1lr
AUCbdvss1
#Pair 6: borderline vs stage II-IV
train$bdvss2 <- mm4_predapp[,2] / (mm4_predapp[,2] +mm4_predapp[,4])
Df = data.frame(p = train$bdvss2, y = train$outcome5, stringsAsFactors = F)
AUCbdvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2) <- c('AUC', 'LL', 'UL')
AUCbdvss2 <- data.frame(AUCbdvss2)
AUCbdvss2lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 4], method = "pepe")
AUCbdvss2[1, 1:3] <- AUCbdvss2lr
AUCbdvss2
#pair 7: borderline vs secondary metastatic
train$bdvssm <- mm4_predapp[,2] / (mm4_predapp[,2] +mm4_predapp[,5])
Df = data.frame(p = train$bdvssm, y = train$outcome5, stringsAsFactors = F)
AUCbdvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssm) <- c('AUC', 'LL', 'UL')
AUCbdvssm <- data.frame(AUCbdvssm)
AUCbdvssmlr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 5], method = "pepe")
AUCbdvssm[1, 1:3] <- AUCbdvssmlr
AUCbdvssm
#pair 8: stage I vs stage II-IV
train$s1vss2 <- mm4_predapp[,3] / (mm4_predapp[,3] +mm4_predapp[,4])
Df = data.frame(p = train$s1vss2, y = train$outcome5, stringsAsFactors = F)
AUCs1vss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2) <- c('AUC', 'LL', 'UL')
AUCs1vss2 <- data.frame(AUCs1vss2)
AUCs1vss2lr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 4], method = "pepe")
AUCs1vss2[1, 1:3] <- AUCs1vss2lr
AUCs1vss2
#pair 9: stage I vs secondary metastatic
train$s1vssm <- mm4_predapp[,3] / (mm4_predapp[,3] +mm4_predapp[,5])
Df = data.frame(p = train$s1vssm, y = train$outcome5, stringsAsFactors = F)
AUCs1vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssm) <- c('AUC', 'LL', 'UL')
AUCs1vssm <- data.frame(AUCs1vssm)
AUCs1vssmlr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 5], method = "pepe")
AUCs1vssm[1, 1:3] <- AUCs1vssmlr
AUCs1vssm
#pair 10: stage II-IV vs secondary metastatic
train$s2vssm <- mm4_predapp[,4] / (mm4_predapp[,4] +mm4_predapp[,5])
Df = data.frame(p = train$s2vssm, y = train$outcome5, stringsAsFactors = F)
AUCs2vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssm) <- c('AUC', 'LL', 'UL')
AUCs2vssm <- data.frame(AUCs2vssm)
AUCs2vssmlr<- auc.nonpara.mw(Df$p[Df$y == 4], Df$p[Df$y == 5], method = "pepe")
AUCs2vssm[1, 1:3] <- AUCs2vssmlr
AUCs2vssm


#FIRTH LR
#Pair 1: benign vs borderline
train$bvsb <- mm5_predapp[, 1] / (mm5_predapp[,1] + mm5_predapp[,2])
Df = data.frame(p = train$bvsb, y = train$outcome5, stringsAsFactors = F)
AUCbvsb <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsb) <- c('AUC', 'LL', 'UL')
AUCbvsb <- data.frame(AUCbvsb)
AUCbvsblr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 2], method = "pepe")
AUCbvsb[1, 1:3] <- AUCbvsblr
#Pair 2: benign vs stage I
train$bvss1 <- mm5_predapp[,1] / (mm5_predapp[,1] +mm5_predapp[,3])
Df = data.frame(p = train$bvss1, y = train$outcome5, stringsAsFactors = F)
AUCbvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1) <- c('AUC', 'LL', 'UL')
AUCbvss1 <- data.frame(AUCbvss1)
AUCbvss1lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 3], method = "pepe")
AUCbvss1[1, 1:3] <- AUCbvss1lr
AUCbvss1$logit.AUC <- logit(AUCbvss1$AUC)
AUCbvss1$logit.se  <- (logit(AUCbvss1$AUC) - logit(AUCbvss1$LL))/1.96
AUCbvss1$logit.var <- AUCbvss1$logit.se^2
AUCbvss1
#Pair 3: benign vs stage II-IV
train$bvss2 <- mm5_predapp[,1] / (mm5_predapp[,1] +mm5_predapp[,4])
Df = data.frame(p = train$bvss2, y = train$outcome5, stringsAsFactors = F)
AUCbvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2) <- c('AUC', 'LL', 'UL')
AUCbvss2 <- data.frame(AUCbvss2)
AUCbvss2lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 4], method = "pepe")
AUCbvss2[1, 1:3] <- AUCbvss2lr
AUCbvss2
#Pair 4: benign vs secondary metastatic
train$bvssm <- mm5_predapp[,1] / (mm5_predapp[,1] +mm5_predapp[,5])
Df = data.frame(p = train$bvssm, y = train$outcome5, stringsAsFactors = F)
AUCbvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssm) <- c('AUC', 'LL', 'UL')
AUCbvssm <- data.frame(AUCbvssm)
AUCbvssmlr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 5], method = "pepe")
AUCbvssm[1, 1:3] <- AUCbvssmlr
AUCbvssm
#Pair 5: borderline vs stage I
train$bdvss1 <- mm5_predapp[,2] / (mm5_predapp[,2] +mm5_predapp[,3])
Df = data.frame(p = train$bdvss1, y = train$outcome5, stringsAsFactors = F)
AUCbdvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1) <- c('AUC', 'LL', 'UL')
AUCbdvss1 <- data.frame(AUCbdvss1)
AUCbdvss1lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 3], method = "pepe")
AUCbdvss1[1, 1:3] <- AUCbdvss1lr
AUCbdvss1
#Pair 6: borderline vs stage II-IV
train$bdvss2 <- mm5_predapp[,2] / (mm5_predapp[,2] +mm5_predapp[,4])
Df = data.frame(p = train$bdvss2, y = train$outcome5, stringsAsFactors = F)
AUCbdvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2) <- c('AUC', 'LL', 'UL')
AUCbdvss2 <- data.frame(AUCbdvss2)
AUCbdvss2lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 4], method = "pepe")
AUCbdvss2[1, 1:3] <- AUCbdvss2lr
AUCbdvss2
#pair 7: borderline vs secondary metastatic
train$bdvssm <- mm5_predapp[,2] / (mm5_predapp[,2] +mm5_predapp[,5])
Df = data.frame(p = train$bdvssm, y = train$outcome5, stringsAsFactors = F)
AUCbdvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssm) <- c('AUC', 'LL', 'UL')
AUCbdvssm <- data.frame(AUCbdvssm)
AUCbdvssmlr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 5], method = "pepe")
AUCbdvssm[1, 1:3] <- AUCbdvssmlr
AUCbdvssm
#pair 8: stage I vs stage II-IV
train$s1vss2 <- mm5_predapp[,3] / (mm5_predapp[,3] +mm5_predapp[,4])
Df = data.frame(p = train$s1vss2, y = train$outcome5, stringsAsFactors = F)
AUCs1vss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2) <- c('AUC', 'LL', 'UL')
AUCs1vss2 <- data.frame(AUCs1vss2)
AUCs1vss2lr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 4], method = "pepe")
AUCs1vss2[1, 1:3] <- AUCs1vss2lr
AUCs1vss2
#pair 9: stage I vs secondary metastatic
train$s1vssm <- mm5_predapp[,3] / (mm5_predapp[,3] +mm5_predapp[,5])
Df = data.frame(p = train$s1vssm, y = train$outcome5, stringsAsFactors = F)
AUCs1vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssm) <- c('AUC', 'LL', 'UL')
AUCs1vssm <- data.frame(AUCs1vssm)
AUCs1vssmlr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 5], method = "pepe")
AUCs1vssm[1, 1:3] <- AUCs1vssmlr
AUCs1vssm
#pair 10: stage II-IV vs secondary metastatic
train$s2vssm <- mm5_predapp[,4] / (mm5_predapp[,4] +mm5_predapp[,5])
Df = data.frame(p = train$s2vssm, y = train$outcome5, stringsAsFactors = F)
AUCs2vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssm) <- c('AUC', 'LL', 'UL')
AUCs2vssm <- data.frame(AUCs2vssm)
AUCs2vssmlr<- auc.nonpara.mw(Df$p[Df$y == 4], Df$p[Df$y == 5], method = "pepe")
AUCs2vssm[1, 1:3] <- AUCs2vssmlr
AUCs2vssm


#RF
#Pair 1: benign vs borderline
train$bvsb <- mm7_predapp[, 1] / (mm7_predapp[,1] + mm7_predapp[,2])
Df = data.frame(p = train$bvsb, y = train$outcome5, stringsAsFactors = F)
AUCbvsb <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsb) <- c('AUC', 'LL', 'UL')
AUCbvsb <- data.frame(AUCbvsb)
AUCbvsblr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 2], method = "pepe")
AUCbvsb[1, 1:3] <- AUCbvsblr
#Pair 2: benign vs stage I
train$bvss1 <- mm7_predapp[,1] / (mm7_predapp[,1] +mm7_predapp[,3])
Df = data.frame(p = train$bvss1, y = train$outcome5, stringsAsFactors = F)
AUCbvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1) <- c('AUC', 'LL', 'UL')
AUCbvss1 <- data.frame(AUCbvss1)
AUCbvss1lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 3], method = "pepe")
AUCbvss1[1, 1:3] <- AUCbvss1lr
AUCbvss1$logit.AUC <- logit(AUCbvss1$AUC)
AUCbvss1$logit.se  <- (logit(AUCbvss1$AUC) - logit(AUCbvss1$LL))/1.96
AUCbvss1$logit.var <- AUCbvss1$logit.se^2
AUCbvss1
#Pair 3: benign vs stage II-IV
train$bvss2 <- mm7_predapp[,1] / (mm7_predapp[,1] +mm7_predapp[,4])
Df = data.frame(p = train$bvss2, y = train$outcome5, stringsAsFactors = F)
AUCbvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2) <- c('AUC', 'LL', 'UL')
AUCbvss2 <- data.frame(AUCbvss2)
AUCbvss2lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 4], method = "pepe")
AUCbvss2[1, 1:3] <- AUCbvss2lr
AUCbvss2
#Pair 4: benign vs secondary metastatic
train$bvssm <- mm7_predapp[,1] / (mm7_predapp[,1] +mm7_predapp[,5])
Df = data.frame(p = train$bvssm, y = train$outcome5, stringsAsFactors = F)
AUCbvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssm) <- c('AUC', 'LL', 'UL')
AUCbvssm <- data.frame(AUCbvssm)
AUCbvssmlr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 5], method = "pepe")
AUCbvssm[1, 1:3] <- AUCbvssmlr
AUCbvssm
#Pair 5: borderline vs stage I
train$bdvss1 <- mm7_predapp[,2] / (mm7_predapp[,2] +mm7_predapp[,3])
Df = data.frame(p = train$bdvss1, y = train$outcome5, stringsAsFactors = F)
AUCbdvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1) <- c('AUC', 'LL', 'UL')
AUCbdvss1 <- data.frame(AUCbdvss1)
AUCbdvss1lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 3], method = "pepe")
AUCbdvss1[1, 1:3] <- AUCbdvss1lr
AUCbdvss1
#Pair 6: borderline vs stage II-IV
train$bdvss2 <- mm7_predapp[,2] / (mm7_predapp[,2] +mm7_predapp[,4])
Df = data.frame(p = train$bdvss2, y = train$outcome5, stringsAsFactors = F)
AUCbdvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2) <- c('AUC', 'LL', 'UL')
AUCbdvss2 <- data.frame(AUCbdvss2)
AUCbdvss2lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 4], method = "pepe")
AUCbdvss2[1, 1:3] <- AUCbdvss2lr
AUCbdvss2
#pair 7: borderline vs secondary metastatic
train$bdvssm <- mm7_predapp[,2] / (mm7_predapp[,2] +mm7_predapp[,5])
Df = data.frame(p = train$bdvssm, y = train$outcome5, stringsAsFactors = F)
AUCbdvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssm) <- c('AUC', 'LL', 'UL')
AUCbdvssm <- data.frame(AUCbdvssm)
AUCbdvssmlr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 5], method = "pepe")
AUCbdvssm[1, 1:3] <- AUCbdvssmlr
AUCbdvssm
#pair 8: stage I vs stage II-IV
train$s1vss2 <- mm7_predapp[,3] / (mm7_predapp[,3] +mm7_predapp[,4])
Df = data.frame(p = train$s1vss2, y = train$outcome5, stringsAsFactors = F)
AUCs1vss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2) <- c('AUC', 'LL', 'UL')
AUCs1vss2 <- data.frame(AUCs1vss2)
AUCs1vss2lr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 4], method = "pepe")
AUCs1vss2[1, 1:3] <- AUCs1vss2lr
AUCs1vss2
#pair 9: stage I vs secondary metastatic
train$s1vssm <- mm7_predapp[,3] / (mm7_predapp[,3] +mm7_predapp[,5])
Df = data.frame(p = train$s1vssm, y = train$outcome5, stringsAsFactors = F)
AUCs1vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssm) <- c('AUC', 'LL', 'UL')
AUCs1vssm <- data.frame(AUCs1vssm)
AUCs1vssmlr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 5], method = "pepe")
AUCs1vssm[1, 1:3] <- AUCs1vssmlr
AUCs1vssm
#pair 10: stage II-IV vs secondary metastatic
train$s2vssm <- mm7_predapp[,4] / (mm7_predapp[,4] +mm7_predapp[,5])
Df = data.frame(p = train$s2vssm, y = train$outcome5, stringsAsFactors = F)
AUCs2vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssm) <- c('AUC', 'LL', 'UL')
AUCs2vssm <- data.frame(AUCs2vssm)
AUCs2vssmlr<- auc.nonpara.mw(Df$p[Df$y == 4], Df$p[Df$y == 5], method = "pepe")
AUCs2vssm[1, 1:3] <- AUCs2vssmlr
AUCs2vssm

#XGBOOST
#Pair 1: benign vs borderline
train$bvsb <- mm8_predapp[, 1] / (mm8_predapp[,1] + mm8_predapp[,2])
Df = data.frame(p = train$bvsb, y = train$outcome5, stringsAsFactors = F)
AUCbvsb <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsb) <- c('AUC', 'LL', 'UL')
AUCbvsb <- data.frame(AUCbvsb)
AUCbvsblr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 2], method = "pepe")
AUCbvsb[1, 1:3] <- AUCbvsblr
#Pair 2: benign vs stage I
train$bvss1 <- mm8_predapp[,1] / (mm8_predapp[,1] +mm8_predapp[,3])
Df = data.frame(p = train$bvss1, y = train$outcome5, stringsAsFactors = F)
AUCbvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1) <- c('AUC', 'LL', 'UL')
AUCbvss1 <- data.frame(AUCbvss1)
AUCbvss1lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 3], method = "pepe")
AUCbvss1[1, 1:3] <- AUCbvss1lr
AUCbvss1$logit.AUC <- logit(AUCbvss1$AUC)
AUCbvss1$logit.se  <- (logit(AUCbvss1$AUC) - logit(AUCbvss1$LL))/1.96
AUCbvss1$logit.var <- AUCbvss1$logit.se^2
AUCbvss1
#Pair 3: benign vs stage II-IV
train$bvss2 <- mm8_predapp[,1] / (mm8_predapp[,1] +mm8_predapp[,4])
Df = data.frame(p = train$bvss2, y = train$outcome5, stringsAsFactors = F)
AUCbvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2) <- c('AUC', 'LL', 'UL')
AUCbvss2 <- data.frame(AUCbvss2)
AUCbvss2lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 4], method = "pepe")
AUCbvss2[1, 1:3] <- AUCbvss2lr
AUCbvss2
#Pair 4: benign vs secondary metastatic
train$bvssm <- mm8_predapp[,1] / (mm8_predapp[,1] +mm8_predapp[,5])
Df = data.frame(p = train$bvssm, y = train$outcome5, stringsAsFactors = F)
AUCbvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssm) <- c('AUC', 'LL', 'UL')
AUCbvssm <- data.frame(AUCbvssm)
AUCbvssmlr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 5], method = "pepe")
AUCbvssm[1, 1:3] <- AUCbvssmlr
AUCbvssm
#Pair 5: borderline vs stage I
train$bdvss1 <- mm8_predapp[,2] / (mm8_predapp[,2] +mm8_predapp[,3])
Df = data.frame(p = train$bdvss1, y = train$outcome5, stringsAsFactors = F)
AUCbdvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1) <- c('AUC', 'LL', 'UL')
AUCbdvss1 <- data.frame(AUCbdvss1)
AUCbdvss1lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 3], method = "pepe")
AUCbdvss1[1, 1:3] <- AUCbdvss1lr
AUCbdvss1
#Pair 6: borderline vs stage II-IV
train$bdvss2 <- mm8_predapp[,2] / (mm8_predapp[,2] +mm8_predapp[,4])
Df = data.frame(p = train$bdvss2, y = train$outcome5, stringsAsFactors = F)
AUCbdvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2) <- c('AUC', 'LL', 'UL')
AUCbdvss2 <- data.frame(AUCbdvss2)
AUCbdvss2lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 4], method = "pepe")
AUCbdvss2[1, 1:3] <- AUCbdvss2lr
AUCbdvss2
#pair 7: borderline vs secondary metastatic
train$bdvssm <- mm8_predapp[,2] / (mm8_predapp[,2] +mm8_predapp[,5])
Df = data.frame(p = train$bdvssm, y = train$outcome5, stringsAsFactors = F)
AUCbdvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssm) <- c('AUC', 'LL', 'UL')
AUCbdvssm <- data.frame(AUCbdvssm)
AUCbdvssmlr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 5], method = "pepe")
AUCbdvssm[1, 1:3] <- AUCbdvssmlr
AUCbdvssm
#pair 8: stage I vs stage II-IV
train$s1vss2 <- mm8_predapp[,3] / (mm8_predapp[,3] +mm8_predapp[,4])
Df = data.frame(p = train$s1vss2, y = train$outcome5, stringsAsFactors = F)
AUCs1vss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2) <- c('AUC', 'LL', 'UL')
AUCs1vss2 <- data.frame(AUCs1vss2)
AUCs1vss2lr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 4], method = "pepe")
AUCs1vss2[1, 1:3] <- AUCs1vss2lr
AUCs1vss2
#pair 9: stage I vs secondary metastatic
train$s1vssm <- mm8_predapp[,3] / (mm8_predapp[,3] +mm8_predapp[,5])
Df = data.frame(p = train$s1vssm, y = train$outcome5, stringsAsFactors = F)
AUCs1vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssm) <- c('AUC', 'LL', 'UL')
AUCs1vssm <- data.frame(AUCs1vssm)
AUCs1vssmlr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 5], method = "pepe")
AUCs1vssm[1, 1:3] <- AUCs1vssmlr
AUCs1vssm
#pair 10: stage II-IV vs secondary metastatic
train$s2vssm <- mm8_predapp[,4] / (mm8_predapp[,4] +mm8_predapp[,5])
Df = data.frame(p = train$s2vssm, y = train$outcome5, stringsAsFactors = F)
AUCs2vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssm) <- c('AUC', 'LL', 'UL')
AUCs2vssm <- data.frame(AUCs2vssm)
AUCs2vssmlr<- auc.nonpara.mw(Df$p[Df$y == 4], Df$p[Df$y == 5], method = "pepe")
AUCs2vssm[1, 1:3] <- AUCs2vssmlr
AUCs2vssm

## NEURAL NETWORKS
#Pair 1: benign vs borderline
train$bvsb <- mm9_predapp[, 1] / (mm9_predapp[,1] + mm9_predapp[,2])
Df = data.frame(p = train$bvsb, y = train$outcome5, stringsAsFactors = F)
AUCbvsb <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsb) <- c('AUC', 'LL', 'UL')
AUCbvsb <- data.frame(AUCbvsb)
AUCbvsblr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 2], method = "pepe")
AUCbvsb[1, 1:3] <- AUCbvsblr
#Pair 2: benign vs stage I
train$bvss1 <- mm9_predapp[,1] / (mm9_predapp[,1] +mm9_predapp[,3])
Df = data.frame(p = train$bvss1, y = train$outcome5, stringsAsFactors = F)
AUCbvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1) <- c('AUC', 'LL', 'UL')
AUCbvss1 <- data.frame(AUCbvss1)
AUCbvss1lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 3], method = "pepe")
AUCbvss1[1, 1:3] <- AUCbvss1lr
AUCbvss1$logit.AUC <- logit(AUCbvss1$AUC)
AUCbvss1$logit.se  <- (logit(AUCbvss1$AUC) - logit(AUCbvss1$LL))/1.96
AUCbvss1$logit.var <- AUCbvss1$logit.se^2
AUCbvss1
#Pair 3: benign vs stage II-IV
train$bvss2 <- mm9_predapp[,1] / (mm9_predapp[,1] +mm9_predapp[,4])
Df = data.frame(p = train$bvss2, y = train$outcome5, stringsAsFactors = F)
AUCbvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2) <- c('AUC', 'LL', 'UL')
AUCbvss2 <- data.frame(AUCbvss2)
AUCbvss2lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 4], method = "pepe")
AUCbvss2[1, 1:3] <- AUCbvss2lr
AUCbvss2
#Pair 4: benign vs secondary metastatic
train$bvssm <- mm9_predapp[,1] / (mm9_predapp[,1] +mm9_predapp[,5])
Df = data.frame(p = train$bvssm, y = train$outcome5, stringsAsFactors = F)
AUCbvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssm) <- c('AUC', 'LL', 'UL')
AUCbvssm <- data.frame(AUCbvssm)
AUCbvssmlr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 5], method = "pepe")
AUCbvssm[1, 1:3] <- AUCbvssmlr
AUCbvssm
#Pair 5: borderline vs stage I
train$bdvss1 <- mm9_predapp[,2] / (mm9_predapp[,2] +mm9_predapp[,3])
Df = data.frame(p = train$bdvss1, y = train$outcome5, stringsAsFactors = F)
AUCbdvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1) <- c('AUC', 'LL', 'UL')
AUCbdvss1 <- data.frame(AUCbdvss1)
AUCbdvss1lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 3], method = "pepe")
AUCbdvss1[1, 1:3] <- AUCbdvss1lr
AUCbdvss1
#Pair 6: borderline vs stage II-IV
train$bdvss2 <- mm9_predapp[,2] / (mm9_predapp[,2] +mm9_predapp[,4])
Df = data.frame(p = train$bdvss2, y = train$outcome5, stringsAsFactors = F)
AUCbdvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2) <- c('AUC', 'LL', 'UL')
AUCbdvss2 <- data.frame(AUCbdvss2)
AUCbdvss2lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 4], method = "pepe")
AUCbdvss2[1, 1:3] <- AUCbdvss2lr
AUCbdvss2
#pair 7: borderline vs secondary metastatic
train$bdvssm <- mm9_predapp[,2] / (mm9_predapp[,2] +mm9_predapp[,5])
Df = data.frame(p = train$bdvssm, y = train$outcome5, stringsAsFactors = F)
AUCbdvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssm) <- c('AUC', 'LL', 'UL')
AUCbdvssm <- data.frame(AUCbdvssm)
AUCbdvssmlr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 5], method = "pepe")
AUCbdvssm[1, 1:3] <- AUCbdvssmlr
AUCbdvssm
#pair 8: stage I vs stage II-IV
train$s1vss2 <- mm9_predapp[,3] / (mm9_predapp[,3] +mm9_predapp[,4])
Df = data.frame(p = train$s1vss2, y = train$outcome5, stringsAsFactors = F)
AUCs1vss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2) <- c('AUC', 'LL', 'UL')
AUCs1vss2 <- data.frame(AUCs1vss2)
AUCs1vss2lr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 4], method = "pepe")
AUCs1vss2[1, 1:3] <- AUCs1vss2lr
AUCs1vss2
#pair 9: stage I vs secondary metastatic
train$s1vssm <- mm9_predapp[,3] / (mm9_predapp[,3] +mm9_predapp[,5])
Df = data.frame(p = train$s1vssm, y = train$outcome5, stringsAsFactors = F)
AUCs1vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssm) <- c('AUC', 'LL', 'UL')
AUCs1vssm <- data.frame(AUCs1vssm)
AUCs1vssmlr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 5], method = "pepe")
AUCs1vssm[1, 1:3] <- AUCs1vssmlr
AUCs1vssm
#pair 10: stage II-IV vs secondary metastatic
train$s2vssm <- mm9_predapp[,4] / (mm9_predapp[,4] +mm9_predapp[,5])
Df = data.frame(p = train$s2vssm, y = train$outcome5, stringsAsFactors = F)
AUCs2vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssm) <- c('AUC', 'LL', 'UL')
AUCs2vssm <- data.frame(AUCs2vssm)
AUCs2vssmlr<- auc.nonpara.mw(Df$p[Df$y == 4], Df$p[Df$y == 5], method = "pepe")
AUCs2vssm[1, 1:3] <- AUCs2vssmlr
AUCs2vssm

#SVM
#Pair 1: benign vs borderline
train$bvsb <- mm10_predapp[, 1] / (mm10_predapp[,1] + mm10_predapp[,2])
Df = data.frame(p = train$bvsb, y = train$outcome5, stringsAsFactors = F)
AUCbvsb <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsb) <- c('AUC', 'LL', 'UL')
AUCbvsb <- data.frame(AUCbvsb)
AUCbvsblr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 2], method = "pepe")
AUCbvsb[1, 1:3] <- AUCbvsblr
#Pair 2: benign vs stage I
train$bvss1 <- mm10_predapp[,1] / (mm10_predapp[,1] +mm10_predapp[,3])
Df = data.frame(p = train$bvss1, y = train$outcome5, stringsAsFactors = F)
AUCbvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1) <- c('AUC', 'LL', 'UL')
AUCbvss1 <- data.frame(AUCbvss1)
AUCbvss1lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 3], method = "pepe")
AUCbvss1[1, 1:3] <- AUCbvss1lr
#Pair 3: benign vs stage II-IV
train$bvss2 <- mm10_predapp[,1] / (mm10_predapp[,1] +mm10_predapp[,4])
Df = data.frame(p = train$bvss2, y = train$outcome5, stringsAsFactors = F)
AUCbvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2) <- c('AUC', 'LL', 'UL')
AUCbvss2 <- data.frame(AUCbvss2)
AUCbvss2lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 4], method = "pepe")
AUCbvss2[1, 1:3] <- AUCbvss2lr
AUCbvss2
#Pair 4: benign vs secondary metastatic
train$bvssm <- mm10_predapp[,1] / (mm10_predapp[,1] +mm10_predapp[,5])
Df = data.frame(p = train$bvssm, y = train$outcome5, stringsAsFactors = F)
AUCbvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssm) <- c('AUC', 'LL', 'UL')
AUCbvssm <- data.frame(AUCbvssm)
AUCbvssmlr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 5], method = "pepe")
AUCbvssm[1, 1:3] <- AUCbvssmlr
AUCbvssm
#Pair 5: borderline vs stage I
train$bdvss1 <- mm10_predapp[,2] / (mm10_predapp[,2] +mm10_predapp[,3])
Df = data.frame(p = train$bdvss1, y = train$outcome5, stringsAsFactors = F)
AUCbdvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1) <- c('AUC', 'LL', 'UL')
AUCbdvss1 <- data.frame(AUCbdvss1)
AUCbdvss1lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 3], method = "pepe")
AUCbdvss1[1, 1:3] <- AUCbdvss1lr
AUCbdvss1
#Pair 6: borderline vs stage II-IV
train$bdvss2 <- mm10_predapp[,2] / (mm10_predapp[,2] +mm10_predapp[,4])
Df = data.frame(p = train$bdvss2, y = train$outcome5, stringsAsFactors = F)
AUCbdvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2) <- c('AUC', 'LL', 'UL')
AUCbdvss2 <- data.frame(AUCbdvss2)
AUCbdvss2lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 4], method = "pepe")
AUCbdvss2[1, 1:3] <- AUCbdvss2lr
AUCbdvss2
#pair 7: borderline vs secondary metastatic
train$bdvssm <- mm10_predapp[,2] / (mm10_predapp[,2] +mm10_predapp[,5])
Df = data.frame(p = train$bdvssm, y = train$outcome5, stringsAsFactors = F)
AUCbdvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssm) <- c('AUC', 'LL', 'UL')
AUCbdvssm <- data.frame(AUCbdvssm)
AUCbdvssmlr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 5], method = "pepe")
AUCbdvssm[1, 1:3] <- AUCbdvssmlr
AUCbdvssm
#pair 8: stage I vs stage II-IV
train$s1vss2 <- mm10_predapp[,3] / (mm10_predapp[,3] +mm10_predapp[,4])
Df = data.frame(p = train$s1vss2, y = train$outcome5, stringsAsFactors = F)
AUCs1vss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2) <- c('AUC', 'LL', 'UL')
AUCs1vss2 <- data.frame(AUCs1vss2)
AUCs1vss2lr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 4], method = "pepe")
AUCs1vss2[1, 1:3] <- AUCs1vss2lr
AUCs1vss2
#pair 9: stage I vs secondary metastatic
train$s1vssm <- mm10_predapp[,3] / (mm10_predapp[,3] +mm10_predapp[,5])
Df = data.frame(p = train$s1vssm, y = train$outcome5, stringsAsFactors = F)
AUCs1vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssm) <- c('AUC', 'LL', 'UL')
AUCs1vssm <- data.frame(AUCs1vssm)
AUCs1vssmlr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 5], method = "pepe")
AUCs1vssm[1, 1:3] <- AUCs1vssmlr
AUCs1vssm
#pair 10: stage II-IV vs secondary metastatic
train$s2vssm <- mm10_predapp[,4] / (mm10_predapp[,4] +mm10_predapp[,5])
Df = data.frame(p = train$s2vssm, y = train$outcome5, stringsAsFactors = F)
AUCs2vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssm) <- c('AUC', 'LL', 'UL')
AUCs2vssm <- data.frame(AUCs2vssm)
AUCs2vssmlr<- auc.nonpara.mw(Df$p[Df$y == 4], Df$p[Df$y == 5], method = "pepe")
AUCs2vssm[1, 1:3] <- AUCs2vssmlr
AUCs2vssm

#######################
##EXTERNAL VALIDATION##
#######################
#Combine the smallest centres into one
levels(test$Center)[levels(test$Center)=="IUK"] <-"Other"
levels(test$Center)[levels(test$Center)=="MIT"] <-"Other"
levels(test$Center)[levels(test$Center)=="FLI"] <-"Other"
levels(test$Center)[levels(test$Center)=="NUK"] <-"Other"
#############################
##1. AUC BENIGN VS MALIGNANT#
#############################

#LOGISTIC REGRESSION WITH TRANSFORMATIONS
mm1_pred <- predict(mm1, newdata=test, type="probs")
test$mm1_pred <- rowSums (mm1_pred[,2:5])

mm1auc <-AUC.IOTA(mm1_pred, outcome=outcome1, center=Center, data=test)


#FOREST PLOT
mm1auc$Plot[17, "RRauc"] <- "   "
mm1auc$Plot[19, "RRauc"] <- "        (0.84 to 0.95)"
mm1auc$Plot[19, "RRprev"] <- "   "
mm1auc$Plot[, "RRcenter"] <- c("Centre", "", "Athens, Greece", "Rome, Italy", "Milan 1, Italy", "Malmö, Sweden", "Stockholm, Sweden", "Genk, Belgium", "Leuven, Belgium", 
                               "Monza, Italy","Cagliari, Italy ", "Other", "Milan 2, Italy", "Pamplona, Spain", "Trieste, Italy", "Katowice, Poland", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("mlr auc.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(mm1auc$Plot,
           align = c("l", "c", "c"),
           mean = mm1auc$dataPlot$AUC,
           lower = mm1auc$dataPlot$LL,
           upper = mm1auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm1auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "Multinomial LR: AUC per center",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm1auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm1auc$Performance$AUC)

dev.off()

#Linear MLR
mm2_pred <- predict(mm2, newdata=test, type="probs")
#calculating overall risk of malignancy
test$mm2_pred <- rowSums (mm2_pred[,2:5])

mm2auc <-AUC.IOTA(test$mm2_pred, outcome=test$outcome1, center=test$Center, data=test)

#FOREST PLOT
mm2auc$Plot[17, "RRauc"] <- "   "
mm2auc$Plot[19, "RRauc"] <- "        (0.84 to 0.94)"
mm2auc$Plot[19, "RRprev"] <- "   "
mm2auc$Plot[, "RRcenter"] <- c("Centre", "", "Athens, Greece", "Rome, Italy", "Milan 1, Italy", "Malmö, Sweden", "Stockholm, Sweden", "Genk, Belgium", "Leuven, Belgium", 
                               "Monza, Italy","Cagliari, Italy ", "Other", "Milan 2, Italy", "Pamplona, Spain", "Trieste, Italy", "Katowice, Poland", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("linear mlr auc.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(mm2auc$Plot,
           align = c("l", "c", "c"),
           mean = mm2auc$dataPlot$AUC,
           lower = mm2auc$dataPlot$LL,
           upper = mm2auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm2auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "Multinomial LR w/o transformations: AUC per center",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm2auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm2auc$Performance$AUC)

dev.off()


#RIDGE REGRESSION
mm4_pred <- predict(mm4, x_test_std, type="response", s=lambda)
mm4_pred <- data.frame(mm4_pred)
#calculating overall risk of malignancy
test$mm4_pred <- rowSums (mm4_pred[,2:5])

mm4auc <-AUC.IOTA(test$mm4_pred, outcome=test$outcome1, center=test$Center, data=test)


#FOREST PLOT
mm4auc$Plot[17, "RRauc"] <- "   "
mm4auc$Plot[19, "RRauc"] <- "        (0.85 to 0.93)"
mm4auc$Plot[19, "RRprev"] <- "   "
mm4auc$Plot[, "RRcenter"] <- c("Centre", "", "Athens, Greece", "Rome, Italy", "Milan 1, Italy", "Malmö, Sweden", "Stockholm, Sweden", "Genk, Belgium", "Leuven, Belgium", 
                               "Monza, Italy","Cagliari, Italy ", "Other", "Milan 2, Italy", "Pamplona, Spain", "Trieste, Italy", "Katowice, Poland", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("ridge auc.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(mm4auc$Plot,
           align = c("l", "c", "c"),
           mean = mm4auc$dataPlot$AUC,
           lower = mm4auc$dataPlot$LL,
           upper = mm4auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm4auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "Multinomial ridge LR: AUC per center",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm4auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm4auc$Performance$AUC)

dev.off()


#FIRTH REGRESSION
#calculating overall risk of malignancy
test$mm5_pred <- rowSums (mm5_pred[,2:5])

mm5auc <-AUC.IOTA(mm5_pred, outcome=outcome1, center=Center, data=test)


#FOREST PLOT
mm5auc$Plot[17, "RRauc"] <- "   "
mm5auc$Plot[19, "RRauc"] <- "        (0.84 to 0.95)"
mm5auc$Plot[19, "RRprev"] <- "   "
mm5auc$Plot[, "RRcenter"] <- c("Centre", "", "Athens, Greece", "Rome, Italy", "Milan 1, Italy", "Malmö, Sweden", "Stockholm, Sweden", "Genk, Belgium", "Leuven, Belgium", 
                               "Monza, Italy","Cagliari, Italy ", "Other", "Milan 2, Italy", "Pamplona, Spain", "Trieste, Italy", "Katowice, Poland", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("firth auc.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(mm5auc$Plot,
           align = c("l", "c", "c"),
           mean = mm5auc$dataPlot$AUC,
           lower = mm5auc$dataPlot$LL,
           upper = mm5auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm5auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "Multinomial firth LR: AUC per center",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm5auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm5auc$Performance$AUC)

dev.off()


#RF
mm7_pred <- predict(mm7, newdata=test, type="prob")
test$mm7_pred <- rowSums (mm7_pred[,2:5])
#1. AUC benign vs malignant with meta analysis
mm7auc <-AUC.IOTA(mm7_pred, outcome=outcome1, center=Center, data=test)


#FOREST PLOT
mm7auc$Plot[17, "RRauc"] <- "   "
mm7auc$Plot[19, "RRauc"] <- "        (0.85 to 0.96)"
mm7auc$Plot[19, "RRprev"] <- "   "
mm7auc$Plot[, "RRcenter"] <- c("Centre", "", "Athens, Greece", "Rome, Italy", "Milan 1, Italy", "Malmö, Sweden", "Stockholm, Sweden", "Genk, Belgium", "Leuven, Belgium", 
                               "Monza, Italy","Cagliari, Italy ", "Other", "Milan 2, Italy", "Pamplona, Spain", "Trieste, Italy", "Katowice, Poland", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("rf auc.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(mm7auc$Plot,
           align = c("l", "c", "c"),
           mean = mm7auc$dataPlot$AUC,
           lower = mm7auc$dataPlot$LL,
           upper = mm7auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm7auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "Multinomial RF: AUC per center",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm7auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm7auc$Performance$AUC)

dev.off()

#XGBOOST
mm8_pred <-predict(object = mm8,newdata = x_test,type = "prob")
test$mm8_pred <- rowSums (mm8_pred[,2:5])

#1.AUC benign vs malignant per center
mm8auc <-AUC.IOTA(mm8_pred, outcome=outcome1, center=Center, data=test)

#FOREST PLOT
mm8auc$Plot[17, "RRauc"] <- "   "
mm8auc$Plot[19, "RRauc"] <- "        (0.84 to 0.96)"
mm8auc$Plot[19, "RRprev"] <- "   "
mm8auc$Plot[, "RRcenter"] <- c("Centre", "", "Athens, Greece", "Rome, Italy", "Milan 1, Italy", "Malmö, Sweden", "Stockholm, Sweden", "Genk, Belgium", "Leuven, Belgium", 
                               "Monza, Italy","Cagliari, Italy ", "Other", "Milan 2, Italy", "Pamplona, Spain", "Trieste, Italy", "Katowice, Poland", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("xg auc bis.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(mm8auc$Plot,
           align = c("l", "c", "c"),
           mean = mm8auc$dataPlot$AUC,
           lower = mm8auc$dataPlot$LL,
           upper = mm8auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm8auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "Multinomial XGBoost: AUC per center",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm8auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm8auc$Performance$AUC)

dev.off()

#NN

mm9_pred <-predict(object = mm9,newdata = d_test,type = "prob")
#calculation of overall risk of malignancy
test$mm9_pred <- rowSums (mm9_pred[,2:5])

#1. AUC benign VS malignant, per center then meta analysis
mm9auc <- AUC.IOTA(pred=mm9_pred, outcome=outcome1, center=Center, data=test)


#FOREST PLOT
mm9auc$Plot[17, "RRauc"] <- "   "
mm9auc$Plot[19, "RRauc"] <- "        (0.86 to 0.95)"
mm9auc$Plot[19, "RRprev"] <- "   "
mm9auc$Plot[, "RRcenter"] <- c("Centre", "", "Athens, Greece", "Rome, Italy", "Milan 1, Italy", "Malmö, Sweden", "Stockholm, Sweden", "Genk, Belgium", "Leuven, Belgium", 
                               "Monza, Italy","Cagliari, Italy ", "Other", "Milan 2, Italy", "Pamplona, Spain", "Trieste, Italy", "Katowice, Poland", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("nn auc.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(mm9auc$Plot,
           align = c("l", "c", "c"),
           mean = mm9auc$dataPlot$AUC,
           lower = mm9auc$dataPlot$LL,
           upper = mm9auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm9auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "Multinomial NN: AUC per center",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm9auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm9auc$Performance$AUC)

dev.off()


#SVM
mm10_pred <- predict(mm10, newdata = d_test %>% select(-V1) , type = "prob")
test$mm10_pred <- rowSums (mm10_pred[,2:5])

mm10auc <- AUC.IOTA(pred=mm10_pred, outcome=outcome1, center=Center, data=test)


#FOREST PLOT
mm10auc$Plot[17, "RRauc"] <- "   "
mm10auc$Plot[19, "RRauc"] <- "        (0.86 to 0.91)"
mm10auc$Plot[19, "RRprev"] <- "   "
mm10auc$Plot[, "RRcenter"] <- c("Centre", "", "Athens, Greece", "Rome, Italy", "Milan 1, Italy", "Malmö, Sweden", "Stockholm, Sweden", "Genk, Belgium", "Leuven, Belgium", 
                                "Monza, Italy","Cagliari, Italy ", "Other", "Milan 2, Italy", "Pamplona, Spain", "Trieste, Italy", "Katowice, Poland", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")

tiff("auc svm.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(mm10auc$Plot,
           align = c("l", "c", "c"),
           mean = mm10auc$dataPlot$AUC,
           lower = mm10auc$dataPlot$LL,
           upper = mm10auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm10auc$IncludedCenters)), TRUE, TRUE, TRUE),
           title = "Multinomial SVM: AUC per center",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm10auc$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm10auc$Performance$AUC)

dev.off()


## OVERALL FOREST PLOT PER MODEL
NA.forest <- mm1auc$Performance[1,]
NA.forest <- NA
Summary.AUC <- rbind(NA.forest, mm1auc$Performance[1,], mm4auc$Performance[1,], mm5auc$Performance[1,], mm2auc$Performance[1,], mm7auc$Performance[1,], mm8auc$Performance[1,], mm9auc$Performance[1,], mm10auc$Performance[1,])
Summary.AUC$Model <- c('', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'RF', "XGBoost", "NN", "SVM")
Summary.AUC.PI <- rbind(NA.forest, mm1auc$Performance[2,], mm4auc$Performance[2,], mm5auc$Performance[2,], mm2auc$Performance[2,], mm7auc$Performance[2,], mm8auc$Performance[2,], mm9auc$Performance[2,], mm10auc$Performance[2,])
Summary.AUC.PI$Model <- c('', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'RF', "XGBoost", "NN", "SVM")
tabletext <- cbind(
  c('Model', 'MLR ', 'Ridge MLR',  'Firth MLR', "Linear MLR", "Random forest", "XGBoost", "Neural network", "Support vector machine"),
  c('AUROC (95% CI)', 
    paste(format(round(mm1auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm1auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm1auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm4auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm4auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm4auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm5auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm5auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm5auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm2auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm2auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm2auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm7auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm7auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm7auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm8auc$Performance$AUC[1], 2), nsmall=2), " (", format(round(mm8auc$Performance$LL[1], 2), nsmall=2), " to ", format(round(mm8auc$Performance$UL[1], 2), nsmall=2), ")", sep= ""),
    paste(format(round(mm9auc$Performance$AUC[1], 2), nsmall=2), " (", format(round(mm9auc$Performance$LL[1], 2), nsmall=2), " to ", format(round(mm9auc$Performance$UL[1], 2), nsmall=2), ")", sep=""),
    paste(format(round(mm10auc$Performance$AUC[1], 2), nsmall=2), " (", format(round(mm10auc$Performance$LL[1], 2), nsmall=2), " to ", format(round(mm10auc$Performance$UL[1], 2), nsmall=2), ")", sep="")),
  c('95% PI', 
    paste0("(", format(round(mm1auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm1auc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm4auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm4auc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm5auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm5auc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm2auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm2auc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm7auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm7auc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm8auc$Performance$LL[2], 2), nsmall=2), " to " , format(round(mm8auc$Performance$UL[2], 2), nsmall=2), ")"),
    paste0("(", format(round(mm9auc$Performance$LL[2], 2), nsmall=2), " to " , format(round(mm9auc$Performance$UL[2], 2), nsmall=2), ")"),
    paste0("(", format(round(mm10auc$Performance$LL[2], 2), nsmall=2), " to " , format(round(mm10auc$Performance$UL[2], 2), nsmall=2), ")")))


tiff("summary plot of binary auroc ext val without ca 2.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE, TRUE,  TRUE, TRUE),
           xlab = "AUROC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.85, 1))
dev.off()

#PDI 

#MLR
mm1pdi <- ests(y=test$outcome5, d=mm1_pred, acc="pdi", level=0.95, method="prob", k=5)
mm1pdi

#LINEAR MLR
mm2pdi <- ests(y=test$outcome5, d=mm2_pred, acc="pdi", level=0.95, method="prob", k=5)
mm2pdi

#RIDGE LOGISTIC REGRESSION
mm4pdi <- ests(y=test$outcome5, d=mm4_pred, acc="pdi", level=0.95, method="prob", k=5)
mm4pdi

#FIRTH LOGISTIC REGRESSION
mm5pdi <- ests(y=test$outcome5, d=mm5_pred, acc="pdi", level=0.95, method="prob", k=5)
mm5pdi

#RANDOM FOREST
mm7pdi <- ests(y=test$outcome5, d=mm7_pred, acc="pdi", level=0.95, method="prob", k=5)
mm7pdi

#XGBOOST
mm8pdi <- ests(y=test$outcome5, d=mm8_pred, acc="pdi", level=0.95, method="prob", k=5)
mm8pdi

#NEURAL NETWORK
mm9pdi <- ests(y=test$outcome5, d=mm9_pred, acc="pdi", level=0.95, method="prob", k=5)
mm9pdi
#SVM
mm10pdi <- ests(y=test$outcome5, d=mm10_pred, acc="pdi", level=0.95, method="prob", k=5)
mm10pdi


#PAIRWISE AUC

#MLR
#Pair 1: benign vs borderline
test$bvsb <- mm1_pred[, 1] / (mm1_pred[,1] + mm1_pred[,2])
Df = data.frame(p = test$bvsb, y = test$outcome5, stringsAsFactors = F)
AUCbvsb <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsb) <- c('AUC', 'LL', 'UL')
AUCbvsb <- data.frame(AUCbvsb)
AUCbvsblr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 2], method = "pepe")
AUCbvsb[1, 1:3] <- AUCbvsblr
#Pair 2: benign vs stage I
test$bvss1 <- mm1_pred[,1] / (mm1_pred[,1] +mm1_pred[,3])
Df = data.frame(p = test$bvss1, y = test$outcome5, stringsAsFactors = F)
AUCbvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1) <- c('AUC', 'LL', 'UL')
AUCbvss1 <- data.frame(AUCbvss1)
AUCbvss1lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 3], method = "pepe")
AUCbvss1[1, 1:3] <- AUCbvss1lr
#Pair 3: benign vs stage II-IV
test$bvss2 <- mm1_pred[,1] / (mm1_pred[,1] +mm1_pred[,4])
Df = data.frame(p = test$bvss2, y = test$outcome5, stringsAsFactors = F)
AUCbvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2) <- c('AUC', 'LL', 'UL')
AUCbvss2 <- data.frame(AUCbvss2)
AUCbvss2lr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 4], method = "pepe")
AUCbvss2[1, 1:3] <- AUCbvss2lr
AUCbvss2
#Pair 4: benign vs secondary metastatic
test$bvssm <- mm1_pred[,1] / (mm1_pred[,1] +mm1_pred[,5])
Df = data.frame(p = test$bvssm, y = test$outcome5, stringsAsFactors = F)
AUCbvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssm) <- c('AUC', 'LL', 'UL')
AUCbvssm <- data.frame(AUCbvssm)
AUCbvssmlr<- auc.nonpara.mw(Df$p[Df$y == 1], Df$p[Df$y == 5], method = "pepe")
AUCbvssm[1, 1:3] <- AUCbvssmlr
AUCbvssm
#Pair 5: borderline vs stage I
test$bdvss1 <- mm1_pred[,2] / (mm1_pred[,2] +mm1_pred[,3])
Df = data.frame(p = test$bdvss1, y = test$outcome5, stringsAsFactors = F)
AUCbdvss1 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1) <- c('AUC', 'LL', 'UL')
AUCbdvss1 <- data.frame(AUCbdvss1)
AUCbdvss1lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 3], method = "pepe")
AUCbdvss1[1, 1:3] <- AUCbdvss1lr
AUCbdvss1
#Pair 6: borderline vs stage II-IV
test$bdvss2 <- mm1_pred[,2] / (mm1_pred[,2] +mm1_pred[,4])
Df = data.frame(p = test$bdvss2, y = test$outcome5, stringsAsFactors = F)
AUCbdvss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2) <- c('AUC', 'LL', 'UL')
AUCbdvss2 <- data.frame(AUCbdvss2)
AUCbdvss2lr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 4], method = "pepe")
AUCbdvss2[1, 1:3] <- AUCbdvss2lr
AUCbdvss2
#pair 7: borderline vs secondary metastatic
test$bdvssm <- mm1_pred[,2] / (mm1_pred[,2] +mm1_pred[,5])
Df = data.frame(p = test$bdvssm, y = test$outcome5, stringsAsFactors = F)
AUCbdvssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssm) <- c('AUC', 'LL', 'UL')
AUCbdvssm <- data.frame(AUCbdvssm)
AUCbdvssmlr<- auc.nonpara.mw(Df$p[Df$y == 2], Df$p[Df$y == 5], method = "pepe")
AUCbdvssm[1, 1:3] <- AUCbdvssmlr
AUCbdvssm
#pair 8: stage I vs stage II-IV
test$s1vss2 <- mm1_pred[,3] / (mm1_pred[,3] +mm1_pred[,4])
Df = data.frame(p = test$s1vss2, y = test$outcome5, stringsAsFactors = F)
AUCs1vss2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2) <- c('AUC', 'LL', 'UL')
AUCs1vss2 <- data.frame(AUCs1vss2)
AUCs1vss2lr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 4], method = "pepe")
AUCs1vss2[1, 1:3] <- AUCs1vss2lr
AUCs1vss2
#pair 9: stage I vs secondary metastatic
test$s1vssm <- mm1_pred[,3] / (mm1_pred[,3] +mm1_pred[,5])
Df = data.frame(p = test$s1vssm, y = test$outcome5, stringsAsFactors = F)
AUCs1vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssm) <- c('AUC', 'LL', 'UL')
AUCs1vssm <- data.frame(AUCs1vssm)
AUCs1vssmlr<- auc.nonpara.mw(Df$p[Df$y == 3], Df$p[Df$y == 5], method = "pepe")
AUCs1vssm[1, 1:3] <- AUCs1vssmlr
AUCs1vssm
#pair 10: stage II-IV vs secondary metastatic
test$s2vssm <- mm1_pred[,4] / (mm1_pred[,4] +mm1_pred[,5])
Df = data.frame(p = test$s2vssm, y = test$outcome5, stringsAsFactors = F)
AUCs2vssm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssm) <- c('AUC', 'LL', 'UL')
AUCs2vssm <- data.frame(AUCs2vssm)
AUCs2vssmlr<- auc.nonpara.mw(Df$p[Df$y == 4], Df$p[Df$y == 5], method = "pepe")
AUCs2vssm[1, 1:3] <- AUCs2vssmlr
AUCs2vssm

#LINEAR MLR
#Pair 1: benign vs borderline
test$bvsblr2 <- mm2_pred[, 1] / (mm2_pred[,1] + mm2_pred[,2])
Dfbvsblr2 = data.frame(p = test$bvsblr2, y = test$outcome5, stringsAsFactors = F)
AUCbvsblr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsblr2) <- c('AUC', 'LL', 'UL')
AUCbvsblr2 <- data.frame(AUCbvsblr2)
bvsbAUClr2<- auc.nonpara.mw(Dfbvsblr2$p[Dfbvsblr2$y == 1], Dfbvsblr2$p[Dfbvsblr2$y == 2], method = "pepe")
AUCbvsblr2[1, 1:3] <- bvsbAUClr2
AUCbvsblr2
#Pair 2: benign vs stage I
test$bvss1lr2 <- mm2_pred[,1] / (mm2_pred[,1] +mm2_pred[,3])
Dfbvss1lr2 = data.frame(p = test$bvss1lr2, y = test$outcome5, stringsAsFactors = F)
AUCbvss1lr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1lr2) <- c('AUC', 'LL', 'UL')
AUCbvss1lr2 <- data.frame(AUCbvss1lr2)
bvss1AUClr2<- auc.nonpara.mw(Dfbvss1lr2$p[Dfbvss1lr2$y == 1], Dfbvss1lr2$p[Dfbvss1lr2$y == 3], method = "pepe")
AUCbvss1lr2[1, 1:3] <- bvss1AUClr2
AUCbvss1lr2
#Pair 3: benign vs stage II-IV
test$bvss2lr2 <- mm2_pred[,1] / (mm2_pred[,1] +mm2_pred[,4])
Dfbvss2lr2 = data.frame(p = test$bvss2lr2, y = test$outcome5, stringsAsFactors = F)
AUCbvss2lr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2lr2) <- c('AUC', 'LL', 'UL')
AUCbvss2lr2 <- data.frame(AUCbvss2lr2)
bvss2AUClr2<- auc.nonpara.mw(Dfbvss2lr2$p[Dfbvss2lr2$y == 1], Dfbvss2lr2$p[Dfbvss2lr2$y == 4], method = "pepe")
AUCbvss2lr2[1, 1:3] <- bvss2AUClr2
AUCbvss2lr2
#Pair 4: benign vs secondary metastatic
test$bvssmlr2 <- mm2_pred[,1] / (mm2_pred[,1] +mm2_pred[,5])
Dfbvssmlr2 = data.frame(p = test$bvssmlr2, y = test$outcome5, stringsAsFactors = F)
AUCbvssmlr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssmlr2) <- c('AUC', 'LL', 'UL')
AUCbvssmlr2 <- data.frame(AUCbvssmlr2)
bvssmAUClr2<- auc.nonpara.mw(Dfbvssmlr2$p[Dfbvssmlr2$y == 1], Dfbvssmlr2$p[Dfbvssmlr2$y == 5], method = "pepe")
AUCbvssmlr2[1, 1:3] <- bvssmAUClr2
AUCbvssmlr2
#Pair 5: borderline vs stage I
test$bdvss1lr2 <- mm2_pred[,2] / (mm2_pred[,2] +mm2_pred[,3])
Dfbdvss1lr2 = data.frame(p = test$bdvss1lr2, y = test$outcome5, stringsAsFactors = F)
AUCbdvss1lr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1lr2) <- c('AUC', 'LL', 'UL')
AUCbdvss1lr2 <- data.frame(AUCbdvss1lr2)
bdvss1AUClr2<- auc.nonpara.mw(Dfbdvss1lr2$p[Dfbdvss1lr2$y == 2], Dfbdvss1lr2$p[Dfbdvss1lr2$y == 3], method = "pepe")
AUCbdvss1lr2[1, 1:3] <- bdvss1AUClr2
AUCbdvss1lr2
#Pair 6: borderline vs stage II-IV
test$bdvss2lr2 <- mm2_pred[,2] / (mm2_pred[,2] +mm2_pred[,4])
Dfbdvss2lr2 = data.frame(p = test$bdvss2lr2, y = test$outcome5, stringsAsFactors = F)
AUCbdvss2lr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2lr2) <- c('AUC', 'LL', 'UL')
AUCbdvss2lr2 <- data.frame(AUCbdvss2lr2)
bdvss2AUClr2<- auc.nonpara.mw(Dfbdvss2lr2$p[Dfbdvss2lr2$y == 2], Dfbdvss2lr2$p[Dfbdvss2lr2$y == 4], method = "pepe")
AUCbdvss2lr2[1, 1:3] <- bdvss2AUClr2
AUCbdvss2lr2
#pair 7: borderline vs secondary metastatic
test$bdvssmlr2 <- mm2_pred[,2] / (mm2_pred[,2] +mm2_pred[,5])
Dfbdvssmlr2 = data.frame(p = test$bdvssmlr2, y = test$outcome5, stringsAsFactors = F)
AUCbdvssmlr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssmlr2) <- c('AUC', 'LL', 'UL')
AUCbdvssmlr2 <- data.frame(AUCbdvssmlr2)
bdvssmAUClr2<- auc.nonpara.mw(Dfbdvssmlr2$p[Dfbdvssmlr2$y == 2], Dfbdvssmlr2$p[Dfbdvssmlr2$y == 5], method = "pepe")
AUCbdvssmlr2[1, 1:3] <- bdvssmAUClr2
AUCbdvssmlr2
#pair 8: stage I vs stage II-IV
test$s1vss2lr2 <- mm2_pred[,3] / (mm2_pred[,3] +mm2_pred[,4])
Dfs1vss2lr2 = data.frame(p = test$s1vss2lr2, y = test$outcome5, stringsAsFactors = F)
AUCs1vss2lr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2lr2) <- c('AUC', 'LL', 'UL')
AUCs1vss2lr2 <- data.frame(AUCs1vss2lr2)
s1vss2AUClr2<- auc.nonpara.mw(Dfs1vss2lr2$p[Dfs1vss2lr2$y == 3], Dfs1vss2lr2$p[Dfs1vss2lr2$y == 4], method = "pepe")
AUCs1vss2lr2[1, 1:3] <- s1vss2AUClr2
AUCs1vss2lr2
#pair 9: stage I vs secondary metastatic
test$s1vssmlr2 <- mm2_pred[,3] / (mm2_pred[,3] +mm2_pred[,5])
Dfs1vssmlr2 = data.frame(p = test$s1vssmlr2, y = test$outcome5, stringsAsFactors = F)
AUCs1vssmlr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssmlr2) <- c('AUC', 'LL', 'UL')
AUCs1vssmlr2 <- data.frame(AUCs1vssmlr2)
s1vssmAUClr2<- auc.nonpara.mw(Dfs1vssmlr2$p[Dfs1vssmlr2$y == 3], Dfs1vssmlr2$p[Dfs1vssmlr2$y == 5], method = "pepe")
AUCs1vssmlr2[1, 1:3] <- s1vss2AUClr2
AUCs1vssmlr2
#pair 10: stage II-IV vs secondary metastatic
test$s2vssmlr2 <- mm2_pred[,4] / (mm2_pred[,4] +mm2_pred[,5])
Dfs2vssmlr2 = data.frame(p = test$s2vssmlr2, y = test$outcome5, stringsAsFactors = F)
AUCs2vssmlr2 <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssmlr2) <- c('AUC', 'LL', 'UL')
AUCs2vssmlr2 <- data.frame(AUCs2vssmlr2)
s2vssmAUClr2<- auc.nonpara.mw(Dfs2vssmlr2$p[Dfs2vssmlr2$y == 4], Dfs2vssmlr2$p[Dfs2vssmlr2$y == 5], method = "pepe")
AUCs2vssmlr2[1, 1:3] <- s2vssmAUClr2
AUCs2vssmlr2

#RIDGE
#Pair 1: benign vs borderline
test$bvsbridge <- mm4_pred[, 1] / (mm4_pred[,1] + mm4_pred[,2])
Dfbvsbridge = data.frame(p = test$bvsbridge, y = test$outcome5, stringsAsFactors = F)
AUCbvsbridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsbridge) <- c('AUC', 'LL', 'UL')
AUCbvsbridge <- data.frame(AUCbvsbridge)
bvsbAUCridge<- auc.nonpara.mw(Dfbvsbridge$p[Dfbvsbridge$y == 1], Dfbvsbridge$p[Dfbvsbridge$y == 2], method = "pepe")
AUCbvsbridge[1, 1:3] <- bvsbAUCridge
AUCbvsbridge
#Pair 2: benign vs stage I
test$bvss1ridge <- mm4_pred[,1] / (mm4_pred[,1] +mm4_pred[,3])
Dfbvss1ridge = data.frame(p = test$bvss1ridge, y = test$outcome5, stringsAsFactors = F)
AUCbvss1ridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1ridge) <- c('AUC', 'LL', 'UL')
AUCbvss1ridge <- data.frame(AUCbvss1ridge)
bvss1AUCridge<- auc.nonpara.mw(Dfbvss1ridge$p[Dfbvss1ridge$y == 1], Dfbvss1ridge$p[Dfbvss1ridge$y == 3], method = "pepe")
AUCbvss1ridge[1, 1:3] <- bvss1AUCridge
AUCbvss1ridge
#Pair 3: benign vs stage II-IV
test$bvss2ridge <- mm4_pred[,1] / (mm4_pred[,1] +mm4_pred[,4])
Dfbvss2ridge = data.frame(p = test$bvss2ridge, y = test$outcome5, stringsAsFactors = F)
AUCbvss2ridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2ridge) <- c('AUC', 'LL', 'UL')
AUCbvss2ridge <- data.frame(AUCbvss2ridge)
bvss2AUCridge<- auc.nonpara.mw(Dfbvss2ridge$p[Dfbvss2ridge$y == 1], Dfbvss2ridge$p[Dfbvss2ridge$y == 4], method = "pepe")
AUCbvss2ridge[1, 1:3] <- bvss2AUCridge
AUCbvss2ridge
#Pair 4: benign vs secondary metastatic
test$bvssmridge <- mm4_pred[,1] / (mm4_pred[,1] +mm4_pred[,5])
Dfbvssmridge = data.frame(p = test$bvssmridge, y = test$outcome5, stringsAsFactors = F)
AUCbvssmridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssmridge) <- c('AUC', 'LL', 'UL')
AUCbvssmridge <- data.frame(AUCbvssmridge)
bvssmAUCridge<- auc.nonpara.mw(Dfbvssmridge$p[Dfbvssmridge$y == 1], Dfbvssmridge$p[Dfbvssmridge$y == 5], method = "pepe")
AUCbvssmridge[1, 1:3] <- bvssmAUCridge
AUCbvssmridge
#Pair 5: borderline vs stage I
test$bdvss1ridge <- mm4_pred[,2] / (mm4_pred[,2] +mm4_pred[,3])
Dfbdvss1ridge = data.frame(p = test$bdvss1ridge, y = test$outcome5, stringsAsFactors = F)
AUCbdvss1ridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1ridge) <- c('AUC', 'LL', 'UL')
AUCbdvss1ridge <- data.frame(AUCbdvss1ridge)
bdvss1AUCridge<- auc.nonpara.mw(Dfbdvss1ridge$p[Dfbdvss1ridge$y == 2], Dfbdvss1ridge$p[Dfbdvss1ridge$y == 3], method = "pepe")
AUCbdvss1ridge[1, 1:3] <- bdvss1AUCridge
AUCbdvss1ridge
#Pair 6: borderline vs stage II-IV
test$bdvss2ridge <- mm4_pred[,2] / (mm4_pred[,2] +mm4_pred[,4])
Dfbdvss2ridge = data.frame(p = test$bdvss2ridge, y = test$outcome5, stringsAsFactors = F)
AUCbdvss2ridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2ridge) <- c('AUC', 'LL', 'UL')
AUCbdvss2ridge <- data.frame(AUCbdvss2ridge)
bdvss2AUCridge<- auc.nonpara.mw(Dfbdvss2ridge$p[Dfbdvss2ridge$y == 2], Dfbdvss2ridge$p[Dfbdvss2ridge$y == 4], method = "pepe")
AUCbdvss2ridge[1, 1:3] <- bdvss2AUCridge
AUCbdvss2ridge
#pair 7: borderline vs secondary metastatic
test$bdvssmridge <- mm4_pred[,2] / (mm4_pred[,2] +mm4_pred[,5])
Dfbdvssmridge = data.frame(p = test$bdvssmridge, y = test$outcome5, stringsAsFactors = F)
AUCbdvssmridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssmridge) <- c('AUC', 'LL', 'UL')
AUCbdvssmridge<- data.frame(AUCbdvssmridge)
bdvssmAUCridge<- auc.nonpara.mw(Dfbdvssmridge$p[Dfbdvssmridge$y == 2], Dfbdvssmridge$p[Dfbdvssmridge$y == 5], method = "pepe")
AUCbdvssmridge[1, 1:3] <- bdvssmAUCridge
AUCbdvssmridge
#pair 8: stage I vs stage II-IV
test$s1vss2ridge <- mm4_pred[,3] / (mm4_pred[,3] +mm4_pred[,4])
Dfs1vss2ridge = data.frame(p = test$s1vss2ridge, y = test$outcome5, stringsAsFactors = F)
AUCs1vss2ridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2ridge) <- c('AUC', 'LL', 'UL')
AUCs1vss2ridge <- data.frame(AUCs1vss2ridge)
s1vss2AUCridge<- auc.nonpara.mw(Dfs1vss2ridge$p[Dfs1vss2ridge$y == 3], Dfs1vss2ridge$p[Dfs1vss2ridge$y == 4], method = "pepe")
AUCs1vss2ridge[1, 1:3] <- s1vss2AUCridge
AUCs1vss2ridge
#pair 9: stage I vs secondary metastatic
test$s1vssmridge <- mm4_pred[,3] / (mm4_pred[,3] +mm4_pred[,5])
Dfs1vssmridge = data.frame(p = test$s1vssmridge, y = test$outcome5, stringsAsFactors = F)
AUCs1vssmridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssmridge) <- c('AUC', 'LL', 'UL')
AUCs1vssmridge <- data.frame(AUCs1vssmridge)
s1vssmAUCridge<- auc.nonpara.mw(Dfs1vssmridge$p[Dfs1vssmridge$y == 3], Dfs1vssmridge$p[Dfs1vssmridge$y == 5], method = "pepe")
AUCs1vssmridge[1, 1:3] <- s1vss2AUCridge
AUCs1vssmridge
#pair 10: stage II-IV vs secondary metastatic
test$s2vssmridge <- mm4_pred[,4] / (mm4_pred[,4] +mm4_pred[,5])
Dfs2vssmridge = data.frame(p = test$s2vssmridge, y = test$outcome5, stringsAsFactors = F)
AUCs2vssmridge <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssmridge) <- c('AUC', 'LL', 'UL')
AUCs2vssmridge <- data.frame(AUCs2vssmridge)
s2vssmAUCridge<- auc.nonpara.mw(Dfs2vssmridge$p[Dfs2vssmridge$y == 4], Dfs2vssmridge$p[Dfs2vssmridge$y == 5], method = "pepe")
AUCs2vssmridge[1, 1:3] <- s2vssmAUCridge
AUCs2vssmridge


#### FIRTH LR
#Pair 1: benign vs borderline
test$bvsbfi <- mm5_pred[, 1] / (mm5_pred[,1] + mm5_pred[,2])
Dfbvsbfi = data.frame(p = test$bvsbfi, y = test$outcome5, stringsAsFactors = F)
AUCbvsbfi <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsbfi) <- c('AUC', 'LL', 'UL')
AUCbvsbfi <- data.frame(AUCbvsbfi)
bvsbAUCfi<- auc.nonpara.mw(Dfbvsbfi$p[Dfbvsbfi$y == 1], Dfbvsbfi$p[Dfbvsbfi$y == 2], method = "pepe")
AUCbvsbfi[1, 1:3] <- bvsbAUCfi
AUCbvsbfi
#Pair 2: benign vs stage I
test$bvss1fi <- mm5_pred[,1] / (mm5_pred[,1] +mm5_pred[,3])
Dfbvss1fi = data.frame(p = test$bvss1fi, y = test$outcome5, stringsAsFactors = F)
AUCbvss1fi <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1fi) <- c('AUC', 'LL', 'UL')
AUCbvss1fi <- data.frame(AUCbvss1fi)
bvss1AUCfi<- auc.nonpara.mw(Dfbvss1fi$p[Dfbvss1fi$y == 1], Dfbvss1fi$p[Dfbvss1fi$y == 3], method = "pepe")
AUCbvss1fi[1, 1:3] <- bvss1AUCfi
AUCbvss1fi
#Pair 3: benign vs stage II-IV
test$bvss2fi <- mm5_pred[,1] / (mm5_pred[,1] +mm5_pred[,4])
Dfbvss2fi = data.frame(p = test$bvss2fi, y = test$outcome5, stringsAsFactors = F)
AUCbvss2fi <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2fi) <- c('AUC', 'LL', 'UL')
AUCbvss2fi <- data.frame(AUCbvss2fi)
bvss2AUCfi<- auc.nonpara.mw(Dfbvss2fi$p[Dfbvss2fi$y == 1], Dfbvss2fi$p[Dfbvss2fi$y == 4], method = "pepe")
AUCbvss2fi[1, 1:3] <- bvss2AUCfi
AUCbvss2fi
#Pair 4: benign vs secondary metastatic
test$bvssmfi <- mm5_pred[,1] / (mm5_pred[,1] +mm5_pred[,5])
Dfbvssmfi = data.frame(p = test$bvssmfi, y = test$outcome5, stringsAsFactors = F)
AUCbvssmfi <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssmfi) <- c('AUC', 'LL', 'UL')
AUCbvssmfi <- data.frame(AUCbvssmfi)
bvssmAUCfi<- auc.nonpara.mw(Dfbvssmfi$p[Dfbvssmfi$y == 1], Dfbvssmfi$p[Dfbvssmfi$y == 5], method = "pepe")
AUCbvssmfi[1, 1:3] <- bvssmAUCfi
AUCbvssmfi
#Pair 5: borderline vs stage I
test$bdvss1fi <- mm5_pred[,2] / (mm5_pred[,2] +mm5_pred[,3])
Dfbdvss1fi = data.frame(p = test$bdvss1fi, y = test$outcome5, stringsAsFactors = F)
AUCbdvss1fi <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1fi) <- c('AUC', 'LL', 'UL')
AUCbdvss1fi <- data.frame(AUCbdvss1fi)
bdvss1AUCfi<- auc.nonpara.mw(Dfbdvss1fi$p[Dfbdvss1fi$y == 2], Dfbdvss1fi$p[Dfbdvss1fi$y == 3], method = "pepe")
AUCbdvss1fi[1, 1:3] <- bdvss1AUCfi
AUCbdvss1fi
#Pair 6: borderline vs stage II-IV
test$bdvss2fi <- mm5_pred[,2] / (mm5_pred[,2] +mm5_pred[,4])
Dfbdvss2fi = data.frame(p = test$bdvss2fi, y = test$outcome5, stringsAsFactors = F)
AUCbdvss2fi <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2fi) <- c('AUC', 'LL', 'UL')
AUCbdvss2fi <- data.frame(AUCbdvss2fi)
bdvss2AUCfi<- auc.nonpara.mw(Dfbdvss2fi$p[Dfbdvss2fi$y == 2], Dfbdvss2fi$p[Dfbdvss2fi$y == 4], method = "pepe")
AUCbdvss2fi[1, 1:3] <- bdvss2AUCfi
AUCbdvss2fi
#pair 7: borderline vs secondary metastatic
test$bdvssmfi <- mm5_pred[,2] / (mm5_pred[,2] +mm5_pred[,5])
Dfbdvssmfi = data.frame(p = test$bdvssmfi, y = test$outcome5, stringsAsFactors = F)
AUCbdvssmfi <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssmfi) <- c('AUC', 'LL', 'UL')
AUCbdvssmfi <- data.frame(AUCbdvssmfi)
bdvssmAUCfi<- auc.nonpara.mw(Dfbdvssmfi$p[Dfbdvssmfi$y == 2], Dfbdvssmfi$p[Dfbdvssmfi$y == 5], method = "pepe")
AUCbdvssmfi[1, 1:3] <- bdvssmAUCfi
AUCbdvssmfi
#pair 8: stage I vs stage II-IV
test$s1vss2fi <- mm5_pred[,3] / (mm5_pred[,3] +mm5_pred[,4])
Dfs1vss2fi = data.frame(p = test$s1vss2fi, y = test$outcome5, stringsAsFactors = F)
AUCs1vss2fi <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2fi) <- c('AUC', 'LL', 'UL')
AUCs1vss2fi <- data.frame(AUCs1vss2fi)
s1vss2AUCfi<- auc.nonpara.mw(Dfs1vss2fi$p[Dfs1vss2fi$y == 3], Dfs1vss2fi$p[Dfs1vss2fi$y == 4], method = "pepe")
AUCs1vss2fi[1, 1:3] <- s1vss2AUCfi
AUCs1vss2fi
#pair 9: stage I vs secondary metastatic
test$s1vssmfi <- mm5_pred[,3] / (mm5_pred[,3] +mm5_pred[,5])
Dfs1vssmfi = data.frame(p = test$s1vssmfi, y = test$outcome5, stringsAsFactors = F)
AUCs1vssmfi <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssmfi) <- c('AUC', 'LL', 'UL')
AUCs1vssmfi <- data.frame(AUCs1vssmfi)
s1vssmAUCfi<- auc.nonpara.mw(Dfs1vssmfi$p[Dfs1vssmfi$y == 3], Dfs1vssmfi$p[Dfs1vssmfi$y == 5], method = "pepe")
AUCs1vssmfi[1, 1:3] <- s1vss2AUCfi
AUCs1vssmfi
#pair 10: stage II-IV vs secondary metastatic
test$s2vssmfi <- mm5_pred[,4] / (mm5_pred[,4] +mm5_pred[,5])
Dfs2vssmfi = data.frame(p = test$s2vssmfi, y = test$outcome5, stringsAsFactors = F)
AUCs2vssmfi <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssmfi) <- c('AUC', 'LL', 'UL')
AUCs2vssmfi <- data.frame(AUCs2vssmfi)
s2vssmAUCfi<- auc.nonpara.mw(Dfs2vssmfi$p[Dfs2vssmfi$y == 4], Dfs2vssmfi$p[Dfs2vssmfi$y == 5], method = "pepe")
AUCs2vssmfi[1, 1:3] <- s2vssmAUCfi
AUCs2vssmfi


#RF
#Pair 1: benign vs borderline
test$bvsbrf <- mm7_pred[, 1] / (mm7_pred[,1] + mm7_pred[,2])
Dfbvsbrf = data.frame(p = test$bvsbrf, y = test$outcome5, stringsAsFactors = F)
AUCbvsbrf <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsbrf) <- c('AUC', 'LL', 'UL')
AUCbvsbrf <- data.frame(AUCbvsbrf)
bvsbAUCrf<- auc.nonpara.mw(Dfbvsbrf$p[Dfbvsbrf$y == 1], Dfbvsbrf$p[Dfbvsbrf$y == 2], method = "pepe")
AUCbvsbrf[1, 1:3] <- bvsbAUCrf
AUCbvsbrf
#Pair 2: benign vs stage I
test$bvss1rf <- mm7_pred[,1] / (mm7_pred[,1] +mm7_pred[,3])
Dfbvss1rf = data.frame(p = test$bvss1rf, y = test$outcome5, stringsAsFactors = F)
AUCbvss1rf <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1rf) <- c('AUC', 'LL', 'UL')
AUCbvss1rf <- data.frame(AUCbvss1rf)
bvss1AUCrf<- auc.nonpara.mw(Dfbvss1rf$p[Dfbvss1rf$y == 1], Dfbvss1rf$p[Dfbvss1rf$y == 3], method = "pepe")
AUCbvss1rf[1, 1:3] <- bvss1AUCrf
AUCbvss1rf
#Pair 3: benign vs stage II-IV
test$bvss2rf <- mm7_pred[,1] / (mm7_pred[,1] +mm7_pred[,4])
Dfbvss2rf = data.frame(p = test$bvss2rf, y = test$outcome5, stringsAsFactors = F)
AUCbvss2rf <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2rf) <- c('AUC', 'LL', 'UL')
AUCbvss2rf <- data.frame(AUCbvss2rf)
bvss2AUCrf<- auc.nonpara.mw(Dfbvss2rf$p[Dfbvss2rf$y == 1], Dfbvss2rf$p[Dfbvss2rf$y == 4], method = "pepe")
AUCbvss2rf[1, 1:3] <- bvss2AUCrf
AUCbvss2rf
#Pair 4: benign vs secondary metastasis
test$bvssmrf <- mm7_pred[,1] / (mm7_pred[,1] +mm7_pred[,5])
Dfbvssmrf = data.frame(p = test$bvssmrf, y = test$outcome5, stringsAsFactors = F)
AUCbvssmrf <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssmrf) <- c('AUC', 'LL', 'UL')
AUCbvssmrf <- data.frame(AUCbvssmrf)
bvssmAUCrf<- auc.nonpara.mw(Dfbvssmrf$p[Dfbvssmrf$y == 1], Dfbvssmrf$p[Dfbvssmrf$y == 5], method = "pepe")
AUCbvssmrf[1, 1:3] <- bvssmAUCrf
AUCbvssmrf
#Pair 5: borderline vs stage I
test$bdvss1rf <- mm7_pred[,2] / (mm7_pred[,2] +mm7_pred[,3])
Dfbdvss1rf = data.frame(p = test$bdvss1rf, y = test$outcome5, stringsAsFactors = F)
AUCbdvss1rf <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1rf) <- c('AUC', 'LL', 'UL')
AUCbdvss1rf <- data.frame(AUCbdvss1rf)
bdvss1AUCrf<- auc.nonpara.mw(Dfbdvss1rf$p[Dfbdvss1rf$y == 2], Dfbdvss1rf$p[Dfbdvss1rf$y == 3], method = "pepe")
AUCbdvss1rf[1, 1:3] <- bdvss1AUCrf
AUCbdvss1rf
#Pair 6: borderline vs stage II-IV
test$bdvss2rf <- mm7_pred[,2] / (mm7_pred[,2] +mm7_pred[,4])
Dfbdvss2rf = data.frame(p = test$bdvss2rf, y = test$outcome5, stringsAsFactors = F)
AUCbdvss2rf <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2rf) <- c('AUC', 'LL', 'UL')
AUCbdvss2rf <- data.frame(AUCbdvss2rf)
bdvss2AUCrf<- auc.nonpara.mw(Dfbdvss2rf$p[Dfbdvss2rf$y == 2], Dfbdvss2rf$p[Dfbdvss2rf$y == 4], method = "pepe")
AUCbdvss2rf[1, 1:3] <- bdvss2AUCrf
AUCbdvss2rf
#pair 7: borderline vs secondary metastasis
test$bdvssmrf <- mm7_pred[,2] / (mm7_pred[,2] +mm7_pred[,5])
Dfbdvssmrf = data.frame(p = test$bdvssmrf, y = test$outcome5, stringsAsFactors = F)
AUCbdvssmrf <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssmrf) <- c('AUC', 'LL', 'UL')
AUCbdvssmrf <- data.frame(AUCbdvssmrf)
bdvssmAUCrf<- auc.nonpara.mw(Dfbdvssmrf$p[Dfbdvssmrf$y == 2], Dfbdvssmrf$p[Dfbdvssmrf$y == 5], method = "pepe")
AUCbdvssmrf[1, 1:3] <- bdvssmAUCrf
AUCbdvssmrf
#pair 8: stage I vs stage II-IV
test$s1vss2rf <- mm7_pred[,3] / (mm7_pred[,3] +mm7_pred[,4])
Dfs1vss2rf = data.frame(p = test$s1vss2rf, y = test$outcome5, stringsAsFactors = F)
AUCs1vss2rf <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2rf) <- c('AUC', 'LL', 'UL')
AUCs1vss2rf <- data.frame(AUCs1vss2rf)
s1vss2AUCrf<- auc.nonpara.mw(Dfs1vss2rf$p[Dfs1vss2rf$y == 3], Dfs1vss2rf$p[Dfs1vss2rf$y == 4], method = "pepe")
AUCs1vss2rf[1, 1:3] <- s1vss2AUCrf
AUCs1vss2rf
#pair 9: stage I vs secondary metastasis
test$s1vssmrf <- mm7_pred[,3] / (mm7_pred[,3] +mm7_pred[,5])
Dfs1vssmrf = data.frame(p = test$s1vssmrf, y = test$outcome5, stringsAsFactors = F)
AUCs1vssmrf <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssmrf) <- c('AUC', 'LL', 'UL')
AUCs1vssmrf <- data.frame(AUCs1vssmrf)
s1vssmAUCrf<- auc.nonpara.mw(Dfs1vssmrf$p[Dfs1vssmrf$y == 3], Dfs1vssmrf$p[Dfs1vssmrf$y == 5], method = "pepe")
AUCs1vssmrf[1, 1:3] <- s1vss2AUCrf
AUCs1vssmrf
#pair 10: stage II-IV vs secondary metastasis
test$s2vssmrf <- mm7_pred[,4] / (mm7_pred[,4] +mm7_pred[,5])
Dfs2vssmrf = data.frame(p = test$s2vssmrf, y = test$outcome5, stringsAsFactors = F)
AUCs2vssmrf <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssmrf) <- c('AUC', 'LL', 'UL')
AUCs2vssmrf <- data.frame(AUCs2vssmrf)
s2vssmAUCrf<- auc.nonpara.mw(Dfs2vssmrf$p[Dfs2vssmrf$y == 4], Dfs2vssmrf$p[Dfs2vssmrf$y == 5], method = "pepe")
AUCs2vssmrf[1, 1:3] <- s2vssmAUCrf
AUCs2vssmrf

#XGBOOST
#Pair 1: benign vs borderline
test$bvsbxg <- mm8_pred[, 1] / (mm8_pred[,1] + mm8_pred[,2])
Dfbvsbxg = data.frame(p = test$bvsbxg, y = test$outcome5, stringsAsFactors = F)
AUCbvsbxg <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsbxg) <- c('AUC', 'LL', 'UL')
AUCbvsbxg <- data.frame(AUCbvsbxg)
bvsbAUCxg<- auc.nonpara.mw(Dfbvsbxg$p[Dfbvsbxg$y == 1], Dfbvsbxg$p[Dfbvsbxg$y == 2], method = "pepe")
AUCbvsbxg[1, 1:3] <- bvsbAUCxg
AUCbvsbxg
#Pair 2: benign vs stage I
test$bvss1xg <- mm8_pred[,1] / (mm8_pred[,1] +mm8_pred[,3])
Dfbvss1xg = data.frame(p = test$bvss1xg, y = test$outcome5, stringsAsFactors = F)
AUCbvss1xg <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1xg) <- c('AUC', 'LL', 'UL')
AUCbvss1xg <- data.frame(AUCbvss1xg)
bvss1AUCxg<- auc.nonpara.mw(Dfbvss1xg$p[Dfbvss1xg$y == 1], Dfbvss1xg$p[Dfbvss1xg$y == 3], method = "pepe")
AUCbvss1xg[1, 1:3] <- bvss1AUCxg
AUCbvss1xg
#Pair 3: benign vs stage II-IV
test$bvss2xg <- mm8_pred[,1] / (mm8_pred[,1] +mm8_pred[,4])
Dfbvss2xg = data.frame(p = test$bvss2xg, y = test$outcome5, stringsAsFactors = F)
AUCbvss2xg <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2xg) <- c('AUC', 'LL', 'UL')
AUCbvss2xg <- data.frame(AUCbvss2xg)
bvss2AUCxg<- auc.nonpara.mw(Dfbvss2xg$p[Dfbvss2xg$y == 1], Dfbvss2xg$p[Dfbvss2xg$y == 4], method = "pepe")
AUCbvss2xg[1, 1:3] <- bvss2AUCxg
AUCbvss2xg
#Pair 4: benign vs secondary metastasis
test$bvssmxg <- mm8_pred[,1] / (mm8_pred[,1] +mm8_pred[,5])
Dfbvssmxg = data.frame(p = test$bvssmxg, y = test$outcome5, stringsAsFactors = F)
AUCbvssmxg <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssmxg) <- c('AUC', 'LL', 'UL')
AUCbvssmxg <- data.frame(AUCbvssmxg)
bvssmAUCxg<- auc.nonpara.mw(Dfbvssmxg$p[Dfbvssmxg$y == 1], Dfbvssmxg$p[Dfbvssmxg$y == 5], method = "pepe")
AUCbvssmxg[1, 1:3] <- bvssmAUCxg
AUCbvssmxg
#Pair 5: borderline vs stage I
test$bdvss1xg <- mm8_pred[,2] / (mm8_pred[,2] +mm8_pred[,3])
Dfbdvss1xg = data.frame(p = test$bdvss1xg, y = test$outcome5, stringsAsFactors = F)
AUCbdvss1xg <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1xg) <- c('AUC', 'LL', 'UL')
AUCbdvss1xg <- data.frame(AUCbdvss1xg)
bdvss1AUCxg<- auc.nonpara.mw(Dfbdvss1xg$p[Dfbdvss1xg$y == 2], Dfbdvss1xg$p[Dfbdvss1xg$y == 3], method = "pepe")
AUCbdvss1xg[1, 1:3] <- bdvss1AUCxg
AUCbdvss1xg
#Pair 6: borderline vs stage II-IV
test$bdvss2xg <- mm8_pred[,2] / (mm8_pred[,2] +mm8_pred[,4])
Dfbdvss2xg = data.frame(p = test$bdvss2xg, y = test$outcome5, stringsAsFactors = F)
AUCbdvss2xg <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2xg) <- c('AUC', 'LL', 'UL')
AUCbdvss2xg <- data.frame(AUCbdvss2xg)
bdvss2AUCxg<- auc.nonpara.mw(Dfbdvss2xg$p[Dfbdvss2xg$y == 2], Dfbdvss2xg$p[Dfbdvss2xg$y == 4], method = "pepe")
AUCbdvss2xg[1, 1:3] <- bdvss2AUCxg
AUCbdvss2xg
#pair 7: borderline vs secondary metastasis
test$bdvssmxg <- mm8_pred[,2] / (mm8_pred[,2] +mm8_pred[,5])
Dfbdvssmxg = data.frame(p = test$bdvssmxg, y = test$outcome5, stringsAsFactors = F)
AUCbdvssmxg <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssmxg) <- c('AUC', 'LL', 'UL')
AUCbdvssmxg <- data.frame(AUCbdvssmxg)
bdvssmAUCxg<- auc.nonpara.mw(Dfbdvssmxg$p[Dfbdvssmxg$y == 2], Dfbdvssmxg$p[Dfbdvssmxg$y == 5], method = "pepe")
AUCbdvssmxg[1, 1:3] <- bdvssmAUCxg
AUCbdvssmxg
#pair 8: stage I vs stage II-IV
test$s1vss2xg <- mm8_pred[,3] / (mm8_pred[,3] +mm8_pred[,4])
Dfs1vss2xg = data.frame(p = test$s1vss2xg, y = test$outcome5, stringsAsFactors = F)
AUCs1vss2xg <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2xg) <- c('AUC', 'LL', 'UL')
AUCs1vss2xg <- data.frame(AUCs1vss2xg)
s1vss2AUCxg<- auc.nonpara.mw(Dfs1vss2xg$p[Dfs1vss2xg$y == 3], Dfs1vss2xg$p[Dfs1vss2xg$y == 4], method = "pepe")
AUCs1vss2xg[1, 1:3] <- s1vss2AUCxg
AUCs1vss2xg
#pair 9: stage I vs secondary metastasis
test$s1vssmxg <- mm8_pred[,3] / (mm8_pred[,3] +mm8_pred[,5])
Dfs1vssmxg = data.frame(p = test$s1vssmxg, y = test$outcome5, stringsAsFactors = F)
AUCs1vssmxg <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssmxg) <- c('AUC', 'LL', 'UL')
AUCs1vssmxg <- data.frame(AUCs1vssmxg)
s1vssmAUCxg<- auc.nonpara.mw(Dfs1vssmxg$p[Dfs1vssmxg$y == 3], Dfs1vssmxg$p[Dfs1vssmxg$y == 5], method = "pepe")
AUCs1vssmxg[1, 1:3] <- s1vss2AUCxg
AUCs1vssmxg
#pair 10: stage II-IV vs secondary metastasis
test$s2vssmxg <- mm8_pred[,4] / (mm8_pred[,4] +mm8_pred[,5])
Dfs2vssmxg = data.frame(p = test$s2vssmxg, y = test$outcome5, stringsAsFactors = F)
AUCs2vssmxg <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssmxg) <- c('AUC', 'LL', 'UL')
AUCs2vssmxg <- data.frame(AUCs2vssmxg)
s2vssmAUCxg<- auc.nonpara.mw(Dfs2vssmxg$p[Dfs2vssmxg$y == 4], Dfs2vssmxg$p[Dfs2vssmxg$y == 5], method = "pepe")
AUCs2vssmxg[1, 1:3] <- s2vssmAUCxg
AUCs2vssmxg

## NEURAL NETWORKS
#Pair 1: benign vs borderline
test$bvsbnn <- mm9_pred[, 1] / (mm9_pred[,1] + mm9_pred[,2])
Dfbvsbnn = data.frame(p = test$bvsbnn, y = test$outcome5, stringsAsFactors = F)
AUCbvsbnn <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsbnn) <- c('AUC', 'LL', 'UL')
AUCbvsbnn <- data.frame(AUCbvsbnn)
bvsbAUCnn<- auc.nonpara.mw(Dfbvsbnn$p[Dfbvsbnn$y == 1], Dfbvsbnn$p[Dfbvsbnn$y == 2], method = "pepe")
AUCbvsbnn[1, 1:3] <- bvsbAUCnn
AUCbvsbnn
#Pair 2: benign vs stage I
test$bvss1nn <- mm9_pred[,1] / (mm9_pred[,1] +mm9_pred[,3])
Dfbvss1nn = data.frame(p = test$bvss1nn, y = test$outcome5, stringsAsFactors = F)
AUCbvss1nn <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1nn) <- c('AUC', 'LL', 'UL')
AUCbvss1nn <- data.frame(AUCbvss1nn)
bvss1AUCnn<- auc.nonpara.mw(Dfbvss1nn$p[Dfbvss1nn$y == 1], Dfbvss1nn$p[Dfbvss1nn$y == 3], method = "pepe")
AUCbvss1nn[1, 1:3] <- bvss1AUCnn
AUCbvss1nn
#Pair 3: benign vs stage II-IV
test$bvss2nn <- mm9_pred[,1] / (mm9_pred[,1] +mm9_pred[,4])
Dfbvss2nn = data.frame(p = test$bvss2nn, y = test$outcome5, stringsAsFactors = F)
AUCbvss2nn <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2nn) <- c('AUC', 'LL', 'UL')
AUCbvss2nn <- data.frame(AUCbvss2nn)
bvss2AUCnn<- auc.nonpara.mw(Dfbvss2nn$p[Dfbvss2nn$y == 1], Dfbvss2nn$p[Dfbvss2nn$y == 4], method = "pepe")
AUCbvss2nn[1, 1:3] <- bvss2AUCnn
AUCbvss2nn
#Pair 4: benign vs secondary metastasis
test$bvssmnn <- mm9_pred[,1] / (mm9_pred[,1] +mm9_pred[,5])
Dfbvssmnn = data.frame(p = test$bvssmnn, y = test$outcome5, stringsAsFactors = F)
AUCbvssmnn <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssmnn) <- c('AUC', 'LL', 'UL')
AUCbvssmnn <- data.frame(AUCbvssmnn)
bvssmAUCnn<- auc.nonpara.mw(Dfbvssmnn$p[Dfbvssmnn$y == 1], Dfbvssmnn$p[Dfbvssmnn$y == 5], method = "pepe")
AUCbvssmnn[1, 1:3] <- bvssmAUCnn
AUCbvssmnn
#Pair 5: borderline vs stage I
test$bdvss1nn <- mm9_pred[,2] / (mm9_pred[,2] +mm9_pred[,3])
Dfbdvss1nn = data.frame(p = test$bdvss1nn, y = test$outcome5, stringsAsFactors = F)
AUCbdvss1nn <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1nn) <- c('AUC', 'LL', 'UL')
AUCbdvss1nn <- data.frame(AUCbdvss1nn)
bdvss1AUCnn<- auc.nonpara.mw(Dfbdvss1nn$p[Dfbdvss1nn$y == 2], Dfbdvss1nn$p[Dfbdvss1nn$y == 3], method = "pepe")
AUCbdvss1nn[1, 1:3] <- bdvss1AUCnn
AUCbdvss1nn
#Pair 6: borderline vs stage II-IV
test$bdvss2nn <- mm9_pred[,2] / (mm9_pred[,2] +mm9_pred[,4])
Dfbdvss2nn = data.frame(p = test$bdvss2nn, y = test$outcome5, stringsAsFactors = F)
AUCbdvss2nn <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2nn) <- c('AUC', 'LL', 'UL')
AUCbdvss2nn <- data.frame(AUCbdvss2nn)
bdvss2AUCnn<- auc.nonpara.mw(Dfbdvss2nn$p[Dfbdvss2nn$y == 2], Dfbdvss2nn$p[Dfbdvss2nn$y == 4], method = "pepe")
AUCbdvss2nn[1, 1:3] <- bdvss2AUCnn
AUCbdvss2nn
#pair 7: borderline vs secondary metastasis
test$bdvssmnn <- mm9_pred[,2] / (mm9_pred[,2] +mm9_pred[,5])
Dfbdvssmnn = data.frame(p = test$bdvssmnn, y = test$outcome5, stringsAsFactors = F)
AUCbdvssmnn <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssmnn) <- c('AUC', 'LL', 'UL')
AUCbdvssmnn <- data.frame(AUCbdvssmnn)
bdvssmAUCnn<- auc.nonpara.mw(Dfbdvssmnn$p[Dfbdvssmnn$y == 2], Dfbdvssmnn$p[Dfbdvssmnn$y == 5], method = "pepe")
AUCbdvssmnn[1, 1:3] <- bdvssmAUCnn
AUCbdvssmnn
#pair 8: stage I vs stage II-IV
test$s1vss2nn <- mm9_pred[,3] / (mm9_pred[,3] +mm9_pred[,4])
Dfs1vss2nn = data.frame(p = test$s1vss2nn, y = test$outcome5, stringsAsFactors = F)
AUCs1vss2nn <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2nn) <- c('AUC', 'LL', 'UL')
AUCs1vss2nn <- data.frame(AUCs1vss2nn)
s1vss2AUCnn<- auc.nonpara.mw(Dfs1vss2nn$p[Dfs1vss2nn$y == 3], Dfs1vss2nn$p[Dfs1vss2nn$y == 4], method = "pepe")
AUCs1vss2nn[1, 1:3] <- s1vss2AUCnn
AUCs1vss2nn
#pair 9: stage I vs secondary metastasis
test$s1vssmnn <- mm9_pred[,3] / (mm9_pred[,3] +mm9_pred[,5])
Dfs1vssmnn = data.frame(p = test$s1vssmnn, y = test$outcome5, stringsAsFactors = F)
AUCs1vssmnn <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssmnn) <- c('AUC', 'LL', 'UL')
AUCs1vssmnn <- data.frame(AUCs1vssmnn)
s1vssmAUCnn<- auc.nonpara.mw(Dfs1vssmnn$p[Dfs1vssmnn$y == 3], Dfs1vssmnn$p[Dfs1vssmnn$y == 5], method = "pepe")
AUCs1vssmnn[1, 1:3] <- s1vssmAUCnn
AUCs1vssmnn
#pair 10: stage II-IV vs secondary metastasis
test$s2vssmnn <- mm9_pred[,4] / (mm9_pred[,4] +mm9_pred[,5])
Dfs2vssmnn = data.frame(p = test$s2vssmnn, y = test$outcome5, stringsAsFactors = F)
AUCs2vssmnn <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssmnn) <- c('AUC', 'LL', 'UL')
AUCs2vssmnn <- data.frame(AUCs2vssmnn)
s2vssmAUCnn<- auc.nonpara.mw(Dfs2vssmnn$p[Dfs2vssmnn$y == 4], Dfs2vssmnn$p[Dfs2vssmnn$y == 5], method = "pepe")
AUCs2vssmnn[1, 1:3] <- s2vssmAUCnn
AUCs2vssmnn

#SVM
#Pair 1: benign vs borderline
test$bvsbsvm <- mm10_pred[, 1] / (mm10_pred[,1] + mm10_pred[,2])
Dfbvsbsvm = data.frame(p = test$bvsbsvm, y = test$outcome5, stringsAsFactors = F)
AUCbvsbsvm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvsbsvm) <- c('AUC', 'LL', 'UL')
AUCbvsbsvm <- data.frame(AUCbvsbsvm)
bvsbAUCsvm<- auc.nonpara.mw(Dfbvsbsvm$p[Dfbvsbsvm$y == 1], Dfbvsbsvm$p[Dfbvsbsvm$y == 2], method = "pepe")
AUCbvsbsvm[1, 1:3] <- bvsbAUCsvm
AUCbvsbsvm
#Pair 2: benign vs stage I
test$bvss1svm <- mm10_pred[,1] / (mm10_pred[,1] +mm10_pred[,3])
Dfbvss1svm = data.frame(p = test$bvss1svm, y = test$outcome5, stringsAsFactors = F)
AUCbvss1svm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss1svm) <- c('AUC', 'LL', 'UL')
AUCbvss1svm <- data.frame(AUCbvss1svm)
bvss1AUCsvm<- auc.nonpara.mw(Dfbvss1svm$p[Dfbvss1svm$y == 1], Dfbvss1svm$p[Dfbvss1svm$y == 3], method = "pepe")
AUCbvss1svm[1, 1:3] <- bvss1AUCsvm
AUCbvss1svm
#Pair 3: benign vs stage II-IV
test$bvss2svm <- mm10_pred[,1] / (mm10_pred[,1] +mm10_pred[,4])
Dfbvss2svm = data.frame(p = test$bvss2, y = test$outcome5, stringsAsFactors = F)
AUCbvss2svm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvss2svm) <- c('AUC', 'LL', 'UL')
AUCbvss2svm <- data.frame(AUCbvss2svm)
bvss2AUCsvm<- auc.nonpara.mw(Dfbvss2svm$p[Dfbvss2svm$y == 1], Dfbvss2svm$p[Dfbvss2svm$y == 4], method = "pepe")
AUCbvss2svm[1, 1:3] <- bvss2AUCsvm
AUCbvss2svm
#Pair 4: benign vs secondary metastasis
test$bvssmsvm <- mm10_pred[,1] / (mm10_pred[,1] +mm10_pred[,5])
Dfbvssmsvm = data.frame(p = test$bvssmsvm, y = test$outcome5, stringsAsFactors = F)
AUCbvssmsvm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbvssmsvm) <- c('AUC', 'LL', 'UL')
AUCbvssmsvm <- data.frame(AUCbvssmsvm)
bvssmAUCsvm<- auc.nonpara.mw(Dfbvssmsvm$p[Dfbvssmsvm$y == 1], Dfbvssmsvm$p[Dfbvssmsvm$y == 5], method = "pepe")
AUCbvssmsvm[1, 1:3] <- bvssmAUCsvm
AUCbvssmsvm
#Pair 5: borderline vs stage I
test$bdvss1svm <- mm10_pred[,2] / (mm10_pred[,2] +mm10_pred[,3])
Dfbdvss1svm = data.frame(p = test$bdvss1svm, y = test$outcome5, stringsAsFactors = F)
AUCbdvss1svm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss1svm) <- c('AUC', 'LL', 'UL')
AUCbdvss1svm <- data.frame(AUCbdvss1svm)
bdvss1AUCsvm<- auc.nonpara.mw(Dfbdvss1svm$p[Dfbdvss1svm$y == 2], Dfbdvss1svm$p[Dfbdvss1svm$y == 3], method = "pepe")
AUCbdvss1svm[1, 1:3] <- bdvss1AUCsvm
AUCbdvss1svm
#Pair 6: borderline vs stage II-IV
test$bdvss2svm <- mm10_pred[,2] / (mm10_pred[,2] +mm10_pred[,4])
Dfbdvss2svm = data.frame(p = test$bdvss2svm, y = test$outcome5, stringsAsFactors = F)
AUCbdvss2svm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvss2svm) <- c('AUC', 'LL', 'UL')
AUCbdvss2svm <- data.frame(AUCbdvss2svm)
bdvss2AUCsvm<- auc.nonpara.mw(Dfbdvss2svm$p[Dfbdvss2svm$y == 2], Dfbdvss2svm$p[Dfbdvss2svm$y == 4], method = "pepe")
AUCbdvss2svm[1, 1:3] <- bdvss2AUCsvm
AUCbdvss2svm
#pair 7: borderline vs secondary metastasis
test$bdvssmsvm <- mm10_pred[,2] / (mm10_pred[,2] +mm10_pred[,5])
Dfbdvssmsvm = data.frame(p = test$bdvssmsvm, y = test$outcome5, stringsAsFactors = F)
AUCbdvssmsvm <- matrix(ncol = 3, nrow = 1)
colnames(AUCbdvssmsvm) <- c('AUC', 'LL', 'UL')
AUCbdvssmsvm <- data.frame(AUCbdvssmsvm)
bdvssmAUCsvm<- auc.nonpara.mw(Dfbdvssmsvm$p[Dfbdvssmsvm$y == 2], Dfbdvssmsvm$p[Dfbdvssmsvm$y == 5], method = "pepe")
AUCbdvssmsvm[1, 1:3] <- bdvssmAUCsvm
AUCbdvssmsvm
#pair 8: stage I vs stage II-IV
test$s1vss2svm <- mm10_pred[,3] / (mm10_pred[,3] +mm10_pred[,4])
Dfs1vss2svm = data.frame(p = test$s1vss2svm, y = test$outcome5, stringsAsFactors = F)
AUCs1vss2svm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vss2svm) <- c('AUC', 'LL', 'UL')
AUCs1vss2svm <- data.frame(AUCs1vss2svm)
s1vss2AUCsvm<- auc.nonpara.mw(Dfs1vss2svm$p[Dfs1vss2svm$y == 3], Dfs1vss2svm$p[Dfs1vss2svm$y == 4], method = "pepe")
AUCs1vss2svm[1, 1:3] <- s1vss2AUCsvm
AUCs1vss2svm
#pair 9: stage I vs secondary metastasis
test$s1vssmsvm <- mm10_pred[,3] / (mm10_pred[,3] +mm10_pred[,5])
Dfs1vssmsvm = data.frame(p = test$s1vssmsvm, y = test$outcome5, stringsAsFactors = F)
AUCs1vssmsvm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs1vssmsvm) <- c('AUC', 'LL', 'UL')
AUCs1vssmsvm <- data.frame(AUCs1vssmsvm)
s1vssmAUCsvm<- auc.nonpara.mw(Dfs1vssmsvm$p[Dfs1vssmsvm$y == 3], Dfs1vssmsvm$p[Dfs1vssmsvm$y == 5], method = "pepe")
AUCs1vssmsvm[1, 1:3] <- s1vssmAUCsvm
AUCs1vssmsvm
#pair 10: stage II-IV vs secondary metastasis
test$s2vssmsvm <- mm10_pred[,4] / (mm10_pred[,4] +mm10_pred[,5])
Dfs2vssmsvm = data.frame(p = test$s2vssmsvm, y = test$outcome5, stringsAsFactors = F)
AUCs2vssmsvm <- matrix(ncol = 3, nrow = 1)
colnames(AUCs2vssmsvm) <- c('AUC', 'LL', 'UL')
AUCs2vssmsvm <- data.frame(AUCs2vssmsvm)
s2vssmAUCsvm<- auc.nonpara.mw(Dfs2vssmsvm$p[Dfs2vssmsvm$y == 4], Dfs2vssmsvm$p[Dfs2vssmsvm$y == 5], method = "pepe")
AUCs2vssmsvm[1, 1:3] <- s2vssmAUCsvm
AUCs2vssmsvm


###########################################
#CALIBRATION SLOPE AND INTERCEPT AND CURVE#
###########################################
####Calibration benign vs malignant

test$Center <- factor(test$Center, labels=c("Athens, Greece", "Milan 1, Italy", "Genk, Belgium", "Leuven, Belgium", "Katowice, Poland", "Malmö, Sweden", 
                                            "Milan 2, Italy", "Monza, Italy", "Other", "Pamplona, Spain", "Rome, Italy", "Cagliari, Italy",
                                            "Stockholm, Sweden", "Trieste, Italy"))

IncludedCenters <- unique(testdf$Center)
#selecting line types, colors and width
lty.centers <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2)
col.centers <- c("#440154FF", "#482677FF", "#404788FF", "#39568CFF", "#2D708EFF", "#238A8DFF", "#1F968BFF", "#20A387FF", "#29AF7FFF","",
                 "#73D055FF", "#95D840FF", "#B8DE29FF", "#FDE725FF")
lwd.centers <- c(1,2,1,2,1,2,1,2,1,2,1,2,1,2)
#MLR
CM1 <-RE.ValProb2(p=mm1_pred, y=outcome1, center=Center,  data=test,MethodCL = "Wald", CalibrLines = "centers")
#Linear MLR
CM2 <-RE.ValProb2(p=mm2_pred, y=outcome1, center=Center,  data=test,  MethodCL = "Wald",CalibrLines = "centers") 

#RIDGE LR
CM4 <-RE.ValProb2(p=mm4_pred, y=outcome1, center=Center, data=test, MethodCL = "Wald", CalibrLines = "centers") 

#FIRTH LR
CM5 <-RE.ValProb2(p=mm5_pred, y=outcome1, center=Center,  data=test,MethodCL = "Wald", CalibrLines = "centers") 

#RF
CM7 <-RE.ValProb2(p=mm7_pred, y=outcome1, center=Center, data=test,  MethodCL = "Wald",CalibrLines = "centers") 

#XGBOOST
CM8 <-RE.ValProb2(p=mm8_pred, y=outcome1, center=Center,  data=test,  MethodCL = "Wald",CalibrLines = "centers") 

#NN
CM9 <-RE.ValProb2(p=mm9_pred, y=outcome1, center=Center,  data=test, MethodCL = "Wald",CalibrLines = "centers") 

#SVM
CM10 <-RE.ValProb2(p=mm10_pred, y=outcome1, center=Center,  data=test, MethodCL = "Wald",CalibrLines = "centers") 



## Predict the outcome for the calibration curve
# MLR
p.M1 = seq(min(test$mm1_pred), max(test$mm1_pred), length = 500)
X = cbind(1, logit(p.M1))
FE = fixef(CM1$FitCalSlope)
OverallCal.M1 = inv.logit(X[order(p.M1), ] %*% FE)
p.M1 = p.M1[order(p.M1)]


#Linear MLR 
p.M2 = seq(min(test$mm2_pred), max(test$mm2_pred), length = 500)
X = cbind(1, logit(p.M2))
FE = fixef(CM2$FitCalSlope)
OverallCal.M2 = inv.logit(X[order(p.M2), ] %*% FE)
p.M2 = p.M2[order(p.M2)]


#RIDGE
p.M4 = seq(min(test$mm4_pred), max(test$mm4_pred), length = 500)
X = cbind(1, logit(p.M4))
FE = fixef(CM4$FitCalSlope)
OverallCal.M4 = inv.logit(X[order(p.M4), ] %*% FE)
p.M4 = p.M4[order(p.M4)]

#FIRTH
p.M5 = seq(min(test$mm5_pred), max(test$mm5_pred), length = 500)
X = cbind(1, logit(p.firth))
FE = fixef(CM5$FitCalSlope)
OverallCal.M5 = inv.logit(X[order(p.M5), ] %*% FE)
p.M5 = p.M5[order(p.M5)]


#RF
p.M7 = seq(min(test$mm7_pred), max(test$mm7_pred), length = 500)
X = cbind(1, logit(p.M7))
FE = fixef(CM1$FitCalSlope)
OverallCal.M7 = inv.logit(X[order(p.M7), ] %*% FE)
p.M7 = p.M7[order(p.M7)]

#XGBOOST
p.M8 = seq(min(test$mm8_pred), max(test$mm8_pred), length = 500)
X = cbind(1, logit(p.M8))
FE = fixef(CM1$FitCalSlope)
OverallCal.M8 = inv.logit(X[order(p.M8), ] %*% FE)
p.M8 = p.M8[order(p.M8)]

#NN
p.M9 = seq(min(test$mm9_pred), max(test$mm9_pred), length = 500)
X = cbind(1, logit(p.M9))
FE = fixef(CM1$FitCalSlope)
OverallCal.M9 = inv.logit(X[order(p.M9), ] %*% FE)
p.M9 = p.M9[order(p.M9)]

#SVM
p.M10 = seq(min(test$mm10_pred), max(test$mm10_pred), length = 500)
X = cbind(1, logit(p.M10))
FE = fixef(CM1$FitCalSlope)
OverallCal.M10 = inv.logit(X[order(p.M10), ] %*% FE)
p.M10 = p.M10[order(p.M10)]



table <- matrix(ncol = 3, nrow = 8)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'Random forest', 'XGBoost', 'Neural network', 'Support vector machine')
table[1, 2:3] <- c(paste0(format(round(CM1$Performance[1,1], 2), nsmall = 2), " (", format(round(CM1$Performance[1,2], 2), nsmall = 2), " to ", format(round(CM1$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(CM1$Performance[2,1], 2), nsmall = 2), " (", format(round(CM1$Performance[2,2], 2), nsmall = 2), " to ", format(round(CM1$Performance[2,3], 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(CM4$Performance[1,1], 2), nsmall = 2), " (", format(round(CM4$Performance[1,2], 2), nsmall = 2), " to ", format(round(CM4$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(CM4$Performance[2,1], 2), nsmall = 2), " (", format(round(CM4$Performance[2,2], 2), nsmall = 2), " to ", format(round(CM4$Performance[2,3], 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(CM5$Performance[1,1], 2), nsmall = 2), " (", format(round(CM5$Performance[1,2], 2), nsmall = 2), " to ", format(round(CM5$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(CM5$Performance[2,1], 2), nsmall = 2), " (", format(round(CM5$Performance[2,2], 2), nsmall = 2), " to ", format(round(CM5$Performance[2,3], 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(CM2$Performance[1,1], 2), nsmall = 2), " (", format(round(CM2$Performance[1,2], 2), nsmall = 2), " to ", format(round(CM2$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(CM2$Performance[2,1], 2), nsmall = 2), " (", format(round(CM2$Performance[2,2], 2), nsmall = 2), " to ", format(round(CM2$Performance[2,3], 2), nsmall = 2), ")"))
table[5, 2:3] <- c(paste0(format(round(CM7$Performance[1,1], 2), nsmall = 2), " (", format(round(CM7$Performance[1,2], 2), nsmall = 2), " to ", format(round(CM7$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(CM7$Performance[2,1], 2), nsmall = 2), " (", format(round(CM7$Performance[2,2], 2), nsmall = 2), " to ", format(round(CM7$Performance[2,3], 2), nsmall = 2), ")"))
table[6, 2:3] <- c(paste0(format(round(CM8$Performance[1,1], 2), nsmall = 2), " (", format(round(CM8$Performance[1,2], 2), nsmall = 2), " to ", format(round(CM8$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(CM8$Performance[2,1], 2), nsmall = 2), " (", format(round(CM8$Performance[2,2], 2), nsmall = 2), " to ", format(round(CM8$Performance[2,3], 2), nsmall = 2), ")"))
table[7, 2:3] <- c(paste0(format(round(CM9$Performance[1,1], 2), nsmall = 2), " (", format(round(CM9$Performance[1,2], 2), nsmall = 2), " to ", format(round(CM9$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(CM9$Performance[2,1], 2), nsmall = 2), " (", format(round(CM9$Performance[2,2], 2), nsmall = 2), " to ", format(round(CM9$Performance[2,3], 2), nsmall = 2), ")"))
table[8, 2:3] <- c(paste0(format(round(CM10$Performance[1,1], 2), nsmall = 2), " (", format(round(CM10$Performance[1,2], 2), nsmall = 2), " to ", format(round(CM10$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(CM10$Performance[2,1], 2), nsmall = 2), " (", format(round(CM10$Performance[2,2], 2), nsmall = 2), " to ", format(round(CM10$Performance[2,3], 2), nsmall = 2), ")"))


## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("multinomial overall calibration table 2.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) # cex.lab = 1.5, cex.axis = 1.5
lines(p.M1, OverallCal.M1, lwd = 2, col = "#440154FF")
lines(p.M4, OverallCal.M4, lwd = 2,lty=3, col = "#404788FF")
lines(p.M5, OverallCal.M5, lwd = 2, lty=5, col = "#33638DFF")
lines(p.M2, OverallCal.M2, lwd=2, lty=3, col="#238A8DFF")
lines(p.M7, OverallCal.M7, lwd=2,lty=2,  col="#20A387FF")
lines(p.M8, OverallCal.M8, lwd=2, lty=3, col="#55C667FF")
lines(p.M9, OverallCal.M9, lwd=2, lty=6, col="#95D840FF")
lines(p.M10, OverallCal.M10, lwd=2, col="#FDE725FF")

legend(x = -0.035, y = 1, legend = c( "Ideal", 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'Random forest', 'XGBoost', 'Neural network', 'Support vector machine'),
       col = c("gray50", "#440154FF", "#404788FF",  "#33638DFF","#238A8DFF",  "#20A387FF", "#55C667FF", "#95D840FF", "#FDE725FF"),
       lty = c(2,1,3,5,3,2,3,6,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1)
addtable2plot(x = 0.25, y = -0.02, table = table, display.colnames= TRUE, cex = 0.6) #cex = 1.2
dev.off()



###############################
##FLEXIBLE MULTINOMIAL CURVES##
###############################


k <- 5
r <- 1
dfr <- 2
outcome <- test$outcome5
cols = c("gray50", "red", "limegreen", "darkblue", "cyan", "magenta")

#MLR
LP1 <- logit(mm1_pred[,-1])
p1 <- mm1_pred
probs1 <- split(p1,col(p1))    
lps1 <- split(LP1,col(LP1))
for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps1[[i]]))}
fitnp1<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)  


#Linear MLR
LP2 <- logit(mm2_pred[,-1])
p2 <- mm2_pred
probs2 <- split(p2,col(p2))    
lps2 <- split(LP2,col(LP2))
for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps2[[i]]))}
fitnp2<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)  


#RIDGE LR
mm4_pred<- data.matrix(mm4_pred)
LP4 <- logit(mm4_pred[,-1])
p4 <- mm4_pred
probs4 <- split(p4,col(p4))    
lps4 <- split(LP4,col(LP4))
for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps4[[i]]))}
fitnp4<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)  


#FIRTH LR
LP5 <- logit(mm5_pred[,-1])
p5 <- mm5_pred
probs5 <- split(p5,col(p5))    
lps5 <- split(LP5,col(LP5))
for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps5[[i]]))}
fitnp5<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)  


#RF
mm7_pred <- data.matrix(mm7_pred)
LP7 <- logit(mm7_pred[,-1])
p7 <- mm7_pred
probs7 <- split(p7,col(p7))    
lps7 <- split(LP7,col(LP7))
for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps7[[i]]))}
fitnp7<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)  


#XGBOOST
mm8_pred <- data.matrix(mm8_pred)
LP8 <- logit(mm8_pred[,-1])
p8 <- mm8_pred
probs8 <- split(p8,col(p8))    
lps8 <- split(LP8,col(LP8))
for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps8[[i]]))}
fitnp8<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)  

#NN
mm9_pred <- data.matrix(mm9_pred)
LP9 <- logit(mm9_pred[,-1])
p9 <- mm9_pred
probs9 <- split(p9,col(p9))    
lps9 <- split(LP9,col(LP9))
for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps9[[i]]))}
fitnp9<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)  

#SVM
mm10_pred <- data.matrix(mm10_pred)
LP10 <- logit(mm10_pred[,-1])
p10 <- mm10_pred
probs10 <- split(p10,col(p10))    
lps10 <- split(LP10,col(LP10))
for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps10[[i]]))}
fitnp10<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)  


#plots

dfmlr <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfmlr) <- c('p1', 'fitted', 'Tumour')
dfmlr <- data.frame(dfmlr)
dfmlr$p1 <- unlist(probs1)
dfmlr$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
fittedmlr <- data.frame(fitted(fitnp1))
dfmlr$fitted <- stack(fittedmlr[1:5])
dfmlr$Tumour <- factor(dfmlr$Tumour, levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))
p1 <-ggplot(dfmlr)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability ")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="MLR")
p1

#ridge
dfridge <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfridge) <- c('p1', 'fitted', 'Tumour')
dfridge <- data.frame(dfridge)
dfridge$p1 <- unlist(probs4)
dfridge$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
fittedridge <- data.frame(fitted(fitnp4))
dfridge$fitted <- stack(fittedridge[1:5])
dfridge$Tumour <- factor(dfridge$Tumour,levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))
p2 <-ggplot(dfridge)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=4, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta"))+geom_text(x=0.1, y=1, label="Ridge MLR")
p2

#firth
dffirth <- matrix(ncol = 3, nrow = 2489*5)
colnames(dffirth) <- c('p1', 'fitted', 'Tumour')
dffirth <- data.frame(dffirth)
dffirth$p1 <- unlist(probs5)
dffirth$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
fittedfirth <- data.frame(fitted(fitnp5))
dffirth$fitted <- stack(fittedfirth[1:5])
dffirth$Tumour <- factor(dffirth$Tumour,levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))
p3 <-ggplot(dffirth)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=4, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+xlim(0,1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability ")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="Firth MLR")
p3


#linear MLR
dflin <- matrix(ncol = 3, nrow = 2489*5)
colnames(dflin) <- c('p1', 'fitted', 'Tumour')
dflin <- data.frame(dflin)
dflin$p1 <- unlist(probs2)
dflin$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
fittedlmr <- data.frame(fitted(fitnp2))
dflin$fitted <- stack(fittedlmr[1:5])
dflin$Tumour <- factor(dflin$Tumour,levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))
#ggplot(dfmlr, aes(x=p1, y=fitted$values, colour=group, group=group))+geom_smooth(method="loess", se=F)+geom_abline(intercept =0 , slope = 1)+xlim(0,1)+ylim(0,1)
p4 <-ggplot(dflin)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability ")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="Linear MLR")
p4

#RF
dfrf <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfrf) <- c('p1', 'fitted', 'Tumour')
dfrf <- data.frame(dfrf)
dfrf$p1 <- unlist(probs7)
dfrf$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
fittedrf <- data.frame(fitted(fitnp7))
dfrf$fitted <- stack(fittedrf[1:5])
dfrf$Tumour <- factor(dfrf$Tumour,levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))
p5 <-ggplot(dfrf)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=4, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability ")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="RF")
p5

#xg
dfxg <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfxg) <- c('p1', 'fitted', 'Tumour')
dfxg <- data.frame(dfxg)
dfxg$p1 <- unlist(probs8)
dfxg$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
fittedxg <- data.frame(fitted(fitnp8))
dfxg$fitted <- stack(fittedxg[1:5])
dfxg$Tumour <- factor(dfxg$Tumour,levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))
p6 <-ggplot(dfxg)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=4, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability ")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="XGBoost")
p6

#NN
dfnn <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfnn) <- c('p1', 'fitted', 'Tumour')
dfnn <- data.frame(dfnn)
dfnn$p1 <- unlist(probs9)
dfnn$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
fittednn <- data.frame(fitted(fitnp9))
dfnn$fitted <- stack(fittednn[1:5])
dfnn$Tumour <- factor(dfnn$Tumour,levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))
p7 <-ggplot(dfnn)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=4, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability ")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="NN")
p7

#SVM
dfsvm <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfsvm) <- c('p1', 'fitted', 'Tumour')
dfsvm <- data.frame(dfsvm)
dfsvm$p1 <- unlist(probs10)
dfsvm$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
fittedsvm <- data.frame(fitted(fitnp10))
dfsvm$fitted <- stack(fittedsvm[1:5])
dfsvm$Tumour <- factor(dfsvm$Tumour,levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))
p8 <-ggplot(dfsvm)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=1.5, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability ")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="SVM")+
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),panel.grid.minor = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),panel.border = element_rect(color = "black",fill = NA))
p8


tiff("multinomial calibration without CA125.tiff", width = 22, height = 31, units = "cm", res = 300)
figure <- ggarrange(p1, p2, p3, p4, p5,p6,p7,p8, ncol=2, nrow=4, common.legend=TRUE, legend="bottom", align="hv")
figure
dev.off()

#net benefit

#MLR
CULR1<- DataWinBugs(pred=test$mm1_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`

#linear MLR
CULR2<- DataWinBugs(pred=test$mm2_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`

#ridge
CULR4<- DataWinBugs(pred=test$mm4_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`

#firth
CULR5<- DataWinBugs(pred=test$mm5_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`

#rf
CULR7<- DataWinBugs(pred=test$mm7_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`

#xg
CULR8 <- DataWinBugs(pred=test$mm8_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`

#nn
CULR9<- DataWinBugs(pred=test$mm9_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`

#svm
CULR10<- DataWinBugs(pred=test$mm10_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`

# plot
data <- read_excel("nb without ca.xlsx")
tiff("nb without ca125 .tiff", width = 16, height = 16, units = "cm", res = 300)

plot(data$threshold,data$`NB treatall`, type = "l", col = "gray50", lty=1, xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.01, 0.35), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(data$threshold, data$`NB LR1`, type='l', col = "#3399CC", lty = 2, lwd = 2)
points(data$threshold,data$`NB RIDGE`, type = 'l', col = "#CCFF66", lty = 3, lwd = 3)
points(data$threshold, data$`NB FIRTH`, type = 'l', col = "#FF9933", lty=4,lwd = 1)
points(data$threshold, data$`NB LR2`, type = 'l', col = "#33FFFF", lty=5,lwd = 2)
points(data$threshold, data$`NB RF`, type = 'l', col = "#332288", lty = 6, lwd = 3)
points(data$threshold, data$`NB XGB`, type = 'l', col = "#00FF66", lty=1,lwd = 1)
points(data$threshold, data$`NB NN`, type = 'l', col = "#FF33CC", lty=2,lwd = 2)
points(data$threshold, data$`NB SVM`, type='l', col = "#882255", lty = 3, lwd = 3)
abline(a = 0, b = 0, lty = 4, lwd = 1, col = "gray50")
legend(x=0.05, y=0.15,legend = c("Treat all","MLR", "Ridge MLR", "Firth MLR", "Linear MLR", "Random forest", "XGBoost", "Neural network", "Support vector machine","Treat none"), ncol = 1, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("gray50", "#3399CC", "#CCFF66",  "#FF9933", "#33FFFF", '#332288', '#00FF66','#FF33CC',"#882255", 'gray50'), 
       lty = c(1,2,3,4,5,6,1,2,3,4), lwd = c(1,2,3,1,2,3,1,2,3,1), cex = 0.8, bty = "n")
dev.off()


############################### WITH CA125 #########################



#1. MLR#


ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  MM1 = multinom(outcome5~oncocenter+Agesd+I(propsolsd)+I(propsolsd^2)+papnrlisd+loc10+Shadows+Ascites+loglesdmax+logca125imp, data = tmp)
  
  
  return(MM1)
}, simplify = F, USE.NAMES = T)

ResultsMM1 = do.call("sapply", args = ArgzCalc)
#Prediction on development data

MM1papp <- data.frame(predict(ResultsMM1,newdata=traindf, type="probs"))

traindf$MM1_pred_be <- rowMeans(subset(MM1papp, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
traindf$MM1_pred_bo <- rowMeans(subset(MM1papp, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
traindf$MM1_pred_bs1 <- rowMeans(subset(MM1papp, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
traindf$MM1_pred_bs2 <- rowMeans(subset(MM1papp, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
traindf$MM1_pred_mt <- rowMeans(subset(MM1papp, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)

traindf$MM1_mal <- rowSums(subset(traindf, select=c(MM1_pred_bo, MM1_pred_bs1, MM1_pred_bs2, MM1_pred_mt)))

MM1papp <- data.frame(traindf$MM1_pred_be, traindf$MM1_pred_bo, traindf$MM1_pred_bs1, traindf$MM1_pred_bs2, traindf$MM1_pred_mt)

#Prediction on validation data
dfpred <- data.frame(predict(ResultsMM1, newdata=testdf, type="probs"))

testdf$MM1_pred_be <- rowMeans(subset(dfpred, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
testdf$MM1_pred_bo <- rowMeans(subset(dfpred, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
testdf$MM1_pred_bs1 <- rowMeans(subset(dfpred, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
testdf$MM1_pred_bs2 <- rowMeans(subset(dfpred, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
testdf$MM1_pred_mt <- rowMeans(subset(dfpred, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)

testdf$MM1_mal <- rowSums(subset(testdf, select=c(MM1_pred_bo, MM1_pred_bs1, MM1_pred_bs2, MM1_pred_mt)))

MM1pred <- data.frame(testdf$MM1_pred_be, testdf$MM1_pred_bo, testdf$MM1_pred_bs1, testdf$MM1_pred_bs2, testdf$MM1_pred_mt)


#histograms of estimated risks

tiff("mlr with ca histogram.tiff", width=14, height=14, units="cm", res=300)
mlr.hist <- ggplot(testdf, aes(x=MM1_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(limits = c(-1000, 3000))
mlr.hist

dev.off()


############
#Linear MLR#
############

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  MM2 = multinom(outcome5~oncocenter+Agesd+propsolsd+papnr+loc10+Shadows+Ascites+lesdmax+ca125imp, data = tmp)
  
  
  return(MM2)
}, simplify = F, USE.NAMES = T)

ResultsMM2 = do.call("sapply", args = ArgzCalc)

#prediction on development data

MM2papp <- data.frame(predict(ResultsMM1,newdata=traindf, type="probs"))

traindf$MM2_pred_be <- rowMeans(subset(MM2papp, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
traindf$MM2_pred_bo <- rowMeans(subset(MM2papp, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
traindf$MM2_pred_bs1 <- rowMeans(subset(MM2papp, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
traindf$MM2_pred_bs2 <- rowMeans(subset(MM2papp, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
traindf$MM2_pred_mt <- rowMeans(subset(MM2papp, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)

traindf$MM2_mal <- rowSums(subset(traindf, select=c(MM2_pred_bo, MM2_pred_bs1, MM2_pred_bs2, MM2_pred_mt)))

MM2papp <- data.frame(traindf$MM2_pred_be, traindf$MM2_pred_bo, traindf$MM2_pred_bs1, traindf$MM2_pred_bs2, traindf$MM2_pred_mt)

#prediction on test data

pMM2 <- data.frame(predict(ResultsMM2, newdata=testdf, type="probs"))

#find a way to do it in a loop
testdf$MM2_pred_be <- rowMeans(subset(pMM2, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
testdf$MM2_pred_bo <- rowMeans(subset(pMM2, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
testdf$MM2_pred_bs1 <- rowMeans(subset(pMM2, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
testdf$MM2_pred_bs2 <- rowMeans(subset(pMM2, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
testdf$MM2_pred_mt <- rowMeans(subset(pMM2, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)

testdf$MM2_mal <- rowSums(subset(testdf, select=c(MM2_pred_bo, MM2_pred_bs1, MM2_pred_bs2, MM2_pred_mt)))

MM2pred <- data.frame(testdf$MM2_pred_be, testdf$MM2_pred_bo, testdf$MM2_pred_bs1, testdf$MM2_pred_bs2, testdf$MM2_pred_mt)

#histogram

tiff("linear mlr with ca histogram.tiff", width=14, height=14, units="cm", res=300)
linearmlr.hist <- ggplot(testdf, aes(x=MM2_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(limits = c(-1000, 3000))
linearmlr.hist

dev.off()

#Ridge logistic regression

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  y_train_std <- data.matrix(traindf[traindf$imp == i,"outcome5"])
  x_train_std <- data.matrix(traindf[traindf$imp == i, c("loglesdmax", "Agesd","oncocenter" ,"propsolsd","papnrlisd", "loc10", "Shadows", "Ascites", "propsolsd2", "logca125imp")])
  #cross-validation to tune lambda
  set.seed(234)
  cv_train_std <- cv.glmnet(x_train_std, y_train_std, nfolds=10, type.measure="deviance", alpha=0, family="multinomial")
  lambda <- cv_train_std$lambda.min
  
  
  #Build the model with optimal lambda value
  MM4 <- glmnet(x_train_std, y_train_std, alpha=0, family="multinomial", lambda=lambda)
  
  return(MM4)
}, simplify = F, USE.NAMES = T)

ResultsMM4 = do.call("sapply", args = ArgzCalc)

#prediction on development data

y_train_std <- data.matrix(traindf[,"outcome5"])
x_train_std <- data.matrix(traindf[, c("loglesdmax", "Agesd","oncocenter" ,"propsolsd","papnrlisd", "loc10", "Shadows", "Ascites", "propsolsd2", "logca125imp")])


MM4papp <- data.frame(predict(ResultsMM4, x_train_std, type="response"))

traindf$MM4_pred_be <- rowMeans(subset(MM4papp, select = c(X1.s0, X1.s0.1,X1.s0.2, X1.s0.3,X1.s0.4,X1.s0.5,X1.s0.6,X1.s0.7,X1.s0.8,X1.s0.9 )), na.rm = TRUE)
traindf$MM4_pred_bo <- rowMeans(subset(MM4papp, select = c(X2.s0, X2.s0.1,X2.s0.2, X2.s0.3,X2.s0.4,X2.s0.5,X2.s0.6,X2.s0.7,X2.s0.8,X2.s0.9 )), na.rm = TRUE)
traindf$MM4_pred_bs1 <- rowMeans(subset(MM4papp, select = c(X3.s0, X3.s0.1,X3.s0.2, X3.s0.3,X3.s0.4,X3.s0.5,X3.s0.6,X3.s0.7,X3.s0.8,X3.s0.9 )), na.rm = TRUE)
traindf$MM4_pred_bs2 <- rowMeans(subset(MM4papp, select = c(X4.s0, X4.s0.1,X4.s0.2, X4.s0.3,X4.s0.4,X4.s0.5,X4.s0.6,X4.s0.7,X4.s0.8,X4.s0.9 )), na.rm = TRUE)
traindf$MM4_pred_mt <- rowMeans(subset(MM4papp, select = c(X5.s0, X5.s0.1,X5.s0.2, X5.s0.3,X5.s0.4,X5.s0.5,X5.s0.6,X5.s0.7,X5.s0.8,X5.s0.9 )), na.rm = TRUE)

traindf$MM4_mal <- rowSums(subset(traindf, select=c(MM4_pred_bo, MM4_pred_bs1, MM4_pred_bs2, MM4_pred_mt)))

MM4papp <- data.frame(traindf$MM4_pred_be, traindf$MM4_pred_bo, traindf$MM4_pred_bs1, traindf$MM4_pred_bs2, traindf$MM4_pred_mt)

#prediction on validation data

y_test_std <- data.matrix(testdf[,"outcome5"])
x_test_std <- data.matrix(testdf[, c("loglesdmax", "Agesd","oncocenter" ,"propsolsd","papnrlisd", "loc10", "Shadows", "Ascites", "propsolsd2", "logca125imp")])

pMM4 <- data.frame(predict(ResultsMM4, x_test_std, type="response"))

testdf$MM4_pred_be <- rowMeans(subset(pMM4, select = c(X1.s0, X1.s0.1,X1.s0.2, X1.s0.3,X1.s0.4,X1.s0.5,X1.s0.6,X1.s0.7,X1.s0.8,X1.s0.9 )), na.rm = TRUE)
testdf$MM4_pred_bo <- rowMeans(subset(pMM4, select = c(X2.s0, X2.s0.1,X2.s0.2, X2.s0.3,X2.s0.4,X2.s0.5,X2.s0.6,X2.s0.7,X2.s0.8,X2.s0.9 )), na.rm = TRUE)
testdf$MM4_pred_bs1 <- rowMeans(subset(pMM4, select = c(X3.s0, X3.s0.1,X3.s0.2, X3.s0.3,X3.s0.4,X3.s0.5,X3.s0.6,X3.s0.7,X3.s0.8,X3.s0.9 )), na.rm = TRUE)
testdf$MM4_pred_bs2 <- rowMeans(subset(pMM4, select = c(X4.s0, X4.s0.1,X4.s0.2, X4.s0.3,X4.s0.4,X4.s0.5,X4.s0.6,X4.s0.7,X4.s0.8,X4.s0.9 )), na.rm = TRUE)
testdf$MM4_pred_mt <- rowMeans(subset(pMM4, select = c(X5.s0, X5.s0.1,X5.s0.2, X5.s0.3,X5.s0.4,X5.s0.5,X5.s0.6,X5.s0.7,X5.s0.8,X5.s0.9 )), na.rm = TRUE)

testdf$MM4_mal <- rowSums(subset(testdf, select=c(MM4_pred_bo, MM4_pred_bs1, MM4_pred_bs2, MM4_pred_mt)))

MM4pred <- data.frame(testdf$MM4_pred_be, testdf$MM4_pred_bo, testdf$MM4_pred_bs1, testdf$MM4_pred_bs2, testdf$MM4_pred_mt)

tiff("ridge mlr with ca histogram.tiff", width=14, height=14, units="cm", res=300)
ridgemlr.hist <- ggplot(testdf, aes(x=MM4_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(limits = c(-1000, 3000))
ridgemlr.hist

dev.off()

#Firth LR 

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  MM5 = brmultinom(outcome5~oncocenter+Agesd+I(propsolsd)+I(propsolsd^2)+papnrlisd+loc10+Shadows+Ascites+loglesdmax + logca125imp, data=tmp,ref=1, type="correction")
  MM5 <- pmlr(outcome5~oncocenter+Agesd+I(propsolsd)+I(propsolsd^2)+papnrlisd+loc10+Shadows+Ascites+loglesdmax, data=tmp,alpha = 0.05, penalized = TRUE, method = "likelihood")
  
  return(MM5)
}, simplify = F, USE.NAMES = T)

ResultsMM5 = do.call("sapply", args = ArgzCalc)

#prediction on development data

MM5papp <- data.frame(predict(ResultsMM5,newdata=traindf, type="probs"))

traindf$MM5_pred_be <- rowMeans(subset(MM5papp, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
traindf$MM5_pred_bo <- rowMeans(subset(MM5papp, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
traindf$MM5_pred_bs1 <- rowMeans(subset(MM5papp, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
traindf$MM5_pred_bs2 <- rowMeans(subset(MM5papp, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
traindf$MM5_pred_mt <- rowMeans(subset(MM5papp, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)


traindf$MM5_mal <- rowSums(subset(traindf, select=c(MM5_pred_bo, MM5_pred_bs1, MM5_pred_bs2, MM5_pred_mt)))

MM5papp <- data.frame(traindf$MM5_pred_be, traindf$MM5_pred_bo, traindf$MM5_pred_bs1, traindf$MM5_pred_bs2, traindf$MM5_pred_mt)

#prediction on test data

pMM5 <- data.frame(predict(ResultsMM5, newdata=testdf, type="probs"))

testdf$MM5_pred_be <- rowMeans(subset(pMM5, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
testdf$MM5_pred_bo <- rowMeans(subset(pMM5, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
testdf$MM5_pred_bs1 <- rowMeans(subset(pMM5, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
testdf$MM5_pred_bs2 <- rowMeans(subset(pMM5, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
testdf$MM5_pred_mt <- rowMeans(subset(pMM5, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)

testdf$MM5_mal <- rowSums(subset(testdf, select=c(MM5_pred_bo, MM5_pred_bs1, MM5_pred_bs2, MM5_pred_mt)))

MM5pred <- data.frame(testdf$MM5_pred_be, testdf$MM5_pred_bo, testdf$MM5_pred_bs1, testdf$MM5_pred_bs2, testdf$MM5_pred_mt)

#histogram
tiff("firth mlr with ca histogram.tiff", width=14, height=14, units="cm", res=300)
firthmlr.hist <- ggplot(testdf, aes(x=MM5_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(limits = c(-1000, 3000))
firthmlr.hist

dev.off()


#RANDOM FOREST
control <- trainControl(method="cv", number=10, classProbs=TRUE, summaryFunction=multiClassSummary, search="random")

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  set.seed(1234)
  MM7 <- train(make.names(outcome5)~oncocenter+Age+propsol+papnr+loc10+Shadows+Ascites+lesdmax+ ca125imp, 
               data = tmp, 
               method = "ranger",
               trControl = control,
               tuneLength = 30,
               metric="logLoss")
  
  return(MM7)
}, simplify = F, USE.NAMES = T)
#getModelInfo('ranger')$ranger$grid
ResultsMM7 = do.call("sapply", args = ArgzCalc)

#prediction on development data
MM7papp <- data.frame(predict(ResultsMM7,newdata=traindf,  type="prob"))


traindf$MM7_pred_be <- rowMeans(subset(MM7papp, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
traindf$MM7_pred_bo <- rowMeans(subset(MM7papp, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
traindf$MM7_pred_bs1 <- rowMeans(subset(MM7papp, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
traindf$MM7_pred_bs2 <- rowMeans(subset(MM7papp, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
traindf$MM7_pred_mt <- rowMeans(subset(MM7papp, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)

traindf$MM7_mal <- rowSums(subset(traindf, select=c(MM7_pred_bo, MM7_pred_bs1, MM7_pred_bs2, MM7_pred_mt)))

MM7papp <- data.frame(traindf$MM7_pred_be, traindf$MM7_pred_bo, traindf$MM7_pred_bs1, traindf$MM7_pred_bs2, traindf$MM7_pred_mt)


#prediction on validation data

pMM7 <- data.frame(predict(ResultsMM73, newdata=testdf, type="prob"))

testdf$MM7_pred_be <- rowMeans(subset(pMM7, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
testdf$MM7_pred_bo <- rowMeans(subset(pMM7, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
testdf$MM7_pred_bs1 <- rowMeans(subset(pMM7, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
testdf$MM7_pred_bs2 <- rowMeans(subset(pMM7, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
testdf$MM7_pred_mt <- rowMeans(subset(pMM7, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)

testdf$MM7_mal <- rowSums(subset(testdf, select=c(MM7_pred_bo, MM7_pred_bs1, MM7_pred_bs2, MM7_pred_mt)))

MM7pred <- data.frame(testdf$MM7_pred_be, testdf$MM7_pred_bo, testdf$MM7_pred_bs1, testdf$MM7_pred_bs2, testdf$MM7_pred_mt)

tiff("rf with ca histogram.tiff", width=14, height=14, units="cm", res=300)
rf.hist <- ggplot(testdf, aes(x=MM7_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(limits = c(-1000, 3000))
rf.hist

dev.off()


#XGBOOST
xgb_trcontrol = trainControl(method = "cv", number = 10, allowParallel = TRUE, 
                             verboseIter = FALSE, returnData = FALSE, summaryFunction=multiClassSummary,   
                             classProbs=TRUE, search="random")

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  train1 <- traindf[traindf$imp == i, c("outcome5", "oncocenter", "lesdmax", "papnr", "loc10", "Shadows", "Ascites", "Age", "propsol", "ca125imp")]
  
  dummies <- dummyVars(~., data=train1[,-1])
  c2 <- predict(dummies, train1[, -1])
  d_training <- as.data.frame(cbind(train1$outcome5, c2))
  
  x_train <- xgb.DMatrix(as.matrix(d_training %>% select (-V1)))
  y_train <- as.factor(make.names(d_training$V1))
  
  
  set.seed(045)
  MM8 = train(x_train, y_train, trControl = xgb_trcontrol, 
              method = "xgbTree", tuneLength=30, metric="logLoss")
  return(MM8)
}, simplify = F, USE.NAMES = T)


ResultsMM8 = do.call("sapply", args = ArgzCalc)


#Prediction on the same data
train1 <- traindf[, c("outcome5", "oncocenter", "lesdmax", "papnr", "loc10", "Shadows", "Ascites", "Age", "propsol", "ca125imp")]

dummies <- dummyVars(~., data=train1[,-1])
c2 <- predict(dummies, train1[, -1])
d_training <- as.data.frame(cbind(train1$outcome5, c2))

x_train <- xgb.DMatrix(as.matrix(d_training %>% select (-V1)))
y_train <- as.factor(make.names(d_training$V1))

MM8papp <- data.frame(predict(ResultsM8, newdata=x_train, type="prob"))


traindf$MM8_pred_be <- rowMeans(subset(MM8papp, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9)), na.rm = TRUE)
traindf$MM8_pred_bo <- rowMeans(subset(MM8papp, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9)), na.rm = TRUE)
traindf$MM8_pred_bs1 <- rowMeans(subset(MM8papp, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9)), na.rm = TRUE)
traindf$MM8_pred_bs2 <- rowMeans(subset(MM8papp, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9)), na.rm = TRUE)
traindf$MM8_pred_mt <- rowMeans(subset(MM8papp, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9)), na.rm = TRUE)

traindf$MM8_mal <- rowSums(subset(traindf, select=c(MM8_pred_bo, MM8_pred_bs1, MM8_pred_bs2, MM8_pred_mt)))

MM8pred <- data.frame(traindf$MM8_pred_be, traindf$MM8_pred_bo, traindf$MM8_pred_bs1, traindf$MM8_pred_bs2, traindf$MM8_pred_mt)


#Prediction on the validation data

test1 <-subset(testdf, select=c("outcome5", "oncocenter", "lesdmax", "papnr", "loc10", "Shadows", "Ascites", "Age", "propsol", "ca125imp"))
dummies <- dummyVars(~., data=test1[, -1])
c2 <- predict(dummies, test1[,-1])
d_test <- as.data.frame(cbind(test1$outcome5, c2))
x_test <- xgb.DMatrix(as.matrix(d_test %>% select(-V1)))
y_test <- as.factor(make.names(d_test$V1))


pMM8 <- data.frame(predict(ResultsM8, newdata=x_test, type="prob"))


testdf$MM8_pred_be <- rowMeans(subset(pMM8, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9)), na.rm = TRUE)
testdf$MM8_pred_bo <- rowMeans(subset(pMM8, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9)), na.rm = TRUE)
testdf$MM8_pred_bs1 <- rowMeans(subset(pMM8, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9)), na.rm = TRUE)
testdf$MM8_pred_bs2 <- rowMeans(subset(pMM8, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9)), na.rm = TRUE)
testdf$MM8_pred_mt <- rowMeans(subset(pMM8, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9)), na.rm = TRUE)

testdf$MM8_mal <- rowSums(subset(testdf, select=c(MM8_pred_bo, MM8_pred_bs1, MM8_pred_bs2, MM8_pred_mt)))

pMM8 <- data.frame(testdf$MM8_pred_be, testdf$MM8_pred_bo, testdf$MM8_pred_bs1, testdf$MM8_pred_bs2, testdf$MM8_pred_mt)

tiff("xg with ca histogram.tiff", width=14, height=14, units="cm", res=300)
xg.hist <- ggplot(testdf, aes(x=MM8_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(limits = c(-1000, 3000))
xg.hist

dev.off()

#NN

ctrl <- trainControl(method="cv",   # 10fold cross validation
                     number=10,         
                     summaryFunction=multiClassSummary,
                     search="random",
                     classProbs=TRUE)

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  train1 <- traindf[traindf$imp == i, c("outcome5", "oncocenter", "lesdmaxsd", "papnr", "loc10", "Shadows", "Ascites", "Agesd", "propsolsd", "ca125imp")]
  dummies <- dummyVars(~., data=train1[,-1])
  c2 <- predict(dummies, train1[, -1])
  d_training <- as.data.frame(cbind(train1$outcome5, c2))
  d_training$V1 <- as.factor(d_training$V1)
  
  
  set.seed(789)
  MM9 <- 
    train(make.names(V1)~. ,
          data = d_training,
          method = "nnet",
          trControl = ctrl,
          trace = FALSE,
          maxit = 1000,
          linout = FALSE,
          tuneLength = 30,
          metric="logLoss") 
  return(MM9)
}, simplify = F, USE.NAMES = T)


ResultsMM9 = do.call("sapply", args = ArgzCalc)

#Prediction on development data
train1 <- traindf[, c("outcome5", "oncocenter", "lesdmaxsd", "papnr", "loc10", "Shadows", "Ascites", "Agesd", "propsolsd", "ca125imp")]
dummies <- dummyVars(~., data=train1[,-1])
c2 <- predict(dummies, train1[, -1])
d_training <- as.data.frame(cbind(train1$outcome5, c2))
d_training$V1 <- as.factor(d_training$V1)

MM9papp <- data.frame(predict(ResultsMM9, newdata=d_training, type="prob"))


traindf$MM9_pred_be <- rowMeans(subset(MM9papp, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
traindf$MM9_pred_bo <- rowMeans(subset(MM9papp, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9)), na.rm = TRUE)
traindf$MM9_pred_bs1 <- rowMeans(subset(MM9papp, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9)), na.rm = TRUE)
traindf$MM9_pred_bs2 <- rowMeans(subset(MM9papp, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
traindf$MM9_pred_mt <- rowMeans(subset(MM9papp, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)

traindf$MM9_mal <- rowSums(subset(traindf, select=c(MM9_pred_bo, MM9_pred_bs1, MM9_pred_bs2, MM9_pred_mt)))

MM9papp <- data.frame(traindf$MM9_pred_be, traindf$MM9_pred_bo, traindf$MM9_pred_bs1, traindf$MM9_pred_bs2, traindf$MM9_pred_mt)

#Prediction on test data

test12 <- subset(testdf, select=c("outcome5", "oncocenter", "lesdmaxsd", "papnr", "loc10", "Shadows", "Ascites", "Agesd", "propsolsd", "ca125imp"))
dummies <- dummyVars(~., data=test12[, -1])
c2 <- predict(dummies, test12[,-1])
d_test <- as.data.frame(cbind(test12$outcome5, c2))
d_test$V1 <- as.factor(d_test$V1)


pMM9 <- data.frame(predict(ResultsMM9, newdata=d_test, type="prob"))

testdf$MM9_pred_be <- rowMeans(subset(pMM9, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
testdf$MM9_pred_bo <- rowMeans(subset(pMM9, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
testdf$MM9_pred_bs1 <- rowMeans(subset(pMM9, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
testdf$MM9_pred_bs2 <- rowMeans(subset(pMM9, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9)), na.rm = TRUE)
testdf$MM9_pred_mt <- rowMeans(subset(pMM9, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9)), na.rm = TRUE)

testdf$MM9_mal <- rowSums(subset(testdf, select=c(MM9_pred_bo, MM9_pred_bs1, MM9_pred_bs2, MM9_pred_mt)))

pMM9 <- data.frame(testdf$MM9_pred_be, testdf$MM9_pred_bo, testdf$MM9_pred_bs1, testdf$MM9_pred_bs2, testdf$MM9_pred_mt)

tiff("nn with ca histogram.tiff", width=14, height=14, units="cm", res=300)
nn.hist <- ggplot(testdf, aes(x=MM9_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(limits = c(-1000, 3000))
nn.hist

dev.off()

#SVM

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  train1 <- traindf[traindf$imp == i, c("outcome5", "oncocenter", "lesdmaxsd", "papnr", "loc10", "Shadows", "Ascites", "Agesd", "propsolsd", "ca125imp")]
  dummies <- dummyVars(~., data=train1[,-1])
  c2 <- predict(dummies, train1[, -1])
  d_training <- as.data.frame(cbind(train1$outcome5, c2))
  d_training$V1 <- as.factor(d_training$V1)
  
  
  set.seed(789)
  MM10 <- train(form=make.names(V1)~.,
                method = "svmRadial",# Radial basis function kernel
                tuneLength = 30,
                metric = "logLoss",
                data=d_training,
                trControl=ctrl) 
  
  return(MM10)
}, simplify = F, USE.NAMES = T)

ResultsMM10 = do.call("sapply", args = ArgzCalc)

#prediction on development data
MM10papp <- data.frame(predict(ResultsMM10, newdata=d_training, type="prob"))

#find a way to do it in a loop
traindf$MM10_pred_be <- rowMeans(subset(MM10papp, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9)), na.rm = TRUE)
traindf$MM10_pred_bo <- rowMeans(subset(MM10papp, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9)), na.rm = TRUE)
traindf$MM10_pred_bs1 <- rowMeans(subset(MM10papp, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9)), na.rm = TRUE)
traindf$MM10_pred_bs2 <- rowMeans(subset(MM10papp, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9)), na.rm = TRUE)
traindf$MM10_pred_mt <- rowMeans(subset(MM10papp, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9)), na.rm = TRUE)

traindf$MM10_mal <- rowSums(subset(traindf, select=c(MM10_pred_bo, MM10_pred_bs1, MM10_pred_bs2, MM10_pred_mt)))

MM10papp <- data.frame(traindf$MM10_pred_be, traindf$MM10_pred_bo, traindf$MM10_pred_bs1, traindf$MM10_pred_bs2, traindf$MM10_pred_mt)

#prediction on test data

pMM10 <- data.frame(predict(ResultsMM10, newdata=d_test, type="prob"))

testdf$MM10_pred_be <- rowMeans(subset(pMM10, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
testdf$MM10_pred_bo <- rowMeans(subset(pMM10, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
testdf$MM10_pred_bs1 <- rowMeans(subset(pMM10, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
testdf$MM10_pred_bs2 <- rowMeans(subset(pMM10, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9)), na.rm = TRUE)
testdf$MM10_pred_mt <- rowMeans(subset(pMM10, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9)), na.rm = TRUE)

testdf$MM10_mal <- rowSums(subset(testdf, select=c(MM10_pred_bo, MM10_pred_bs1, MM10_pred_bs2, MM10_pred_mt)))

MM10pred <- data.frame(testdf$MM10_pred_be, testdf$MM10_pred_bo, testdf$MM10_pred_bs1, testdf$MM10_pred_bs2, testdf$MM10_pred_mt)

tiff("svm with ca histogram.tiff", width=14, height=14, units="cm", res=300)
svm.hist <- ggplot(testdf, aes(x=MM10_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(limits = c(-1000, 3000))
svm.hist

dev.off()


#SCATTERPLOTS

#scatterplots of the probability of a benign tumour# 

#combine all predictions of a benign tumour from each model
benigncombined <- matrix(ncol = 8, nrow = 24890)
colnames(benigncombined) <- c('MLR', 'Ridge', 'Firth', 'Linear', 'RF', 'XGBoost', 'NN', 'SVM')
benigncombined <- data.frame(benigncombined)
benigncombined$MLR <- testdf$MM1_pred_be
benigncombined$Ridge <- testdf$MM4_pred_be
benigncombined$Firth <- testdf$MM5_pred_be
benigncombined$Linear <- testdf$MM2_pred_be
benigncombined$RF <- testdf$MM7_pred_be
benigncombined$XGBoost <- testdf$MM8_pred_be
benigncombined$NN <- testdf$MM9_pred_be
benigncombined$SVM <- testdf$MM10_pred_be


plotfun <- function(x,y,...){
  points(x,y,...)
  abline(0,1, col="red")
}

tiff("scatterplot benign with ca multi models benign not combined.tiff", width=14, height=14, units="cm", res=300)
pairs(benigncombined[,1:8], pch=".", upper.panel=NULL, cex=0.5,  cex.labels=0.95,xlim=c(0,1), ylim=c(0,1)) #lower.panel=plotfun,
dev.off()


#scatterplot of estimated risk of a borderline tumour 
bocombined <- matrix(ncol = 8, nrow = 24890)
colnames(bocombined) <- c('MLR', 'Ridge', 'Firth', 'Linear', 'RF', 'XGBoost', 'NN', 'SVM')
bocombined <- data.frame(bocombined)
bocombined$MLR <- testdf$MM1_pred_bo
bocombined$Ridge <- testdf$MM4_pred_bo
bocombined$Firth <- testdf$MM5_pred_bo
bocombined$Linear <- testdf$MM2_pred_bo
bocombined$RF <- testdf$MM7_pred_bo
bocombined$XGBoost <- testdf$MM8_pred_bo
bocombined$NN <- testdf$MM9_pred_bo
bocombined$SVM <- testdf$MM10_pred_bo


tiff("scatterplot borderline with ca 3.tiff", width=14, height=14, units="cm", res=300)
pairs(bocombined[,1:8], pch=".", upper.panel=NULL, cex=0.5, lower.panel=plotfun, cex.labels=0.95, xlim=c(0,1), ylim=c(0,1))
dev.off()

#Scatterplot of estimated risk of a stage I primary invasive tumour
s1combined <- matrix(ncol = 8, nrow = 24890)
colnames(s1combined) <- c('MLR', 'Ridge', 'Firth', 'Linear', 'RF', 'XGBoost', 'NN', 'SVM')
s1combined <- data.frame(s1combined)
s1combined$MLR <- testdf$MM1_pred_bs1
s1combined$Ridge <- testdf$MM4_pred_bs1
s1combined$Firth <- testdf$MM5_pred_bs1
s1combined$Linear <- testdf$MM2_pred_bs1
s1combined$RF <- testdf$MM7_pred_bs1
s1combined$XGBoost <- testdf$MM8_pred_bs1
s1combined$NN <- testdf$MM9_pred_bs1
s1combined$SVM <- testdf$MM10_pred_bs1

tiff("scatterplot stage1 with ca 3.tiff", width=14, height=14, units="cm", res=300)
pairs(s1combined[,1:8], pch=".", upper.panel=NULL, cex=0.5, lower.panel=plotfun, cex.labels=0.95, xlim=c(0,1), ylim=c(0,1))
dev.off()

#Scatterplot of estimated risk of a stage II-IV primary invasive tumour
s2combined <- matrix(ncol = 8, nrow = 24890)
colnames(s2combined) <- c('MLR', 'Ridge', 'Firth', 'Linear', 'RF', 'XGBoost', 'NN', 'SVM')
s2combined <- data.frame(s2combined)
s2combined$MLR <- testdf$MM1_pred_bs2
s2combined$Ridge <- testdf$MM4_pred_bs2
s2combined$Firth <- testdf$MM5_pred_bs2
s2combined$Linear <- testdf$MM2_pred_bs2
s2combined$RF <- testdf$MM7_pred_bs2
s2combined$XGBoost <- testdf$MM8_pred_bs2
s2combined$NN <- testdf$MM9_pred_bs2
s2combined$SVM <- testdf$MM10_pred_bs2


tiff("scatterplot s2 with ca 3.tiff", width=14, height=14, units="cm", res=300)
pairs(s2combined[,1:8], pch=".", upper.panel=NULL, cex=0.5, lower.panel=plotfun, cex.labels=0.95, xlim=c(0,1), ylim=c(0,1))
dev.off()

#Scatterplots of estimated risk of a secondary metastatic tumour
mtcombined <- matrix(ncol = 8, nrow = 24890)
colnames(mtcombined) <- c('MLR', 'Ridge', 'Firth', 'Linear', 'RF', 'XGBoost', 'NN', 'SVM')
mtcombined <- data.frame(mtcombined)
mtcombined$MLR <- testdf$MM1_pred_mt
mtcombined$Ridge <- testdf$MM4_pred_mt
mtcombined$Firth <- testdf$MM5_pred_mt
mtcombined$Linear <- testdf$MM2_pred_mt
mtcombined$RF <- testdf$MM7_pred_mt
mtcombined$XGBoost <- testdf$MM8_pred_mt
mtcombined$NN <- testdf$MM9_pred_mt
mtcombined$SVM <- testdf$MM10_pred_mt


tiff("scatterplot metastatic with ca 3.tiff", width=14, height=14, units="cm", res=300)
pairs(mtcombined[,1:8], pch=".", upper.panel=NULL, cex=0.5, lower.panel=plotfun, cex.labels=0.95, xlim=c(0,1), ylim=c(0,1))
dev.off()


#histograms for each outcome
testdf$outcome55 <- factor(testdf$outcome5, levels = c("1", "2", "3", "4", "5"), labels=c("Benign", "BO", "St.I", "S.II-IV", "SMeta"))

#benign
p1.1 <-ggplot(testdf, aes(x=MM1_pred_be))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a benign outcome (MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"))+
  scale_y_continuous(limits=c(0,7200), breaks=c(0,2000,4000,6000))
p1.2 <-ggplot(testdf, aes(x=MM4_pred_be))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a benign outcome (Ridge MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,7200), breaks=c(0,2000,4000,6000))

p1.3 <-ggplot(testdf, aes(x=MM5_pred_be))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a benign outcome (Firth MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"))+
  scale_y_continuous(limits=c(0,7200), breaks=c(0,2000,4000,6000))

p1.4 <-ggplot(testdf, aes(x=MM2_pred_be))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+xlab("Estimated probability of a benign outcome (Linear MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels = c("0", "0.5", "1"))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,7200), breaks=c(0,2000,4000,6000))

p1.5 <-ggplot(testdf, aes(x=MM7_pred_be))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a benign outcome (RF)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"))+
  scale_y_continuous(limits=c(0,7200), breaks=c(0,2000,4000,6000))

p1.6 <-ggplot(testdf, aes(x=MM8_pred_be))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+xlab("Estimated probability of a benign outcome (XGBoost)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,7200), breaks=c(0,2000,4000,6000))

p1.7 <-ggplot(testdf, aes(x=MM9_pred_be))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a benign outcome (NN)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"))+
  scale_y_continuous(limits=c(0,7200), breaks=c(0,2000,4000,6000))

p1.8 <-ggplot(testdf, aes(x=MM10_pred_be))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+xlab("Estimated probability of a benign outcome (SVM)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,7200), breaks=c(0,2000,4000,6000))

tiff("histograms benign all models with ca125 2.tiff", width = 26, height = 31, units = "cm", res = 300)
figure <- ggarrange(p1.1, p1.2, p1.3, p1.4, p1.5,p1.6,p1.7,p1.8, ncol=2, nrow=4, common.legend=TRUE, legend="bottom", align="hv")
figure
dev.off()


#Borderline
p2.1 <-ggplot(testdf, aes(x=MM1_pred_bo))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a borderline outcome (MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,10500), breaks=c(0,2000,4000,6000, 8000, 10000))

p2.2 <-ggplot(testdf, aes(x=MM4_pred_bo))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a borderline outcome (Ridge MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,10500), breaks=c(0,2000,4000,6000, 8000, 10000))

p2.3 <-ggplot(testdf, aes(x=MM5_pred_bo))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a borderline outcome (Firth MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,10500), breaks=c(0,2000,4000,6000, 8000, 10000))

p2.4 <-ggplot(testdf, aes(x=MM2_pred_bo))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a borderline outcome (Linear MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,10500), breaks=c(0,2000,4000,6000, 8000, 10000))

p2.5 <-ggplot(testdf, aes(x=MM7_pred_bo))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a borderline outcome (RF)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,10500), breaks=c(0,2000,4000,6000, 8000, 10000))

p2.6 <-ggplot(testdf, aes(x=MM8_pred_bo))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+xlab("Estimated probability of a borderline outcome (XGBoost)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,10500), breaks=c(0,2000,4000,6000, 8000, 10000))

p2.7 <-ggplot(testdf, aes(x=MM9_pred_bo))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a borderline outcome (NN)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,10500), breaks=c(0,2000,4000,6000, 8000, 10000))

p2.8 <-ggplot(testdf, aes(x=MM10_pred_bo))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+xlab("Estimated probability of a borderline outcome (SVM)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,10500), breaks=c(0,2000,4000,6000, 8000, 10000))

tiff("histograms borderline all models with ca125 2.tiff", width = 25, height = 31, units = "cm", res = 300)
figure <- ggarrange(p2.1, p2.2, p2.3, p2.4, p2.5,p2.6,p2.7,p2.8, ncol=2, nrow=4, common.legend=TRUE, legend="bottom", align="hv")
figure
dev.off()

#stage 1
p3.1 <-ggplot(testdf, aes(x=MM1_pred_bs1))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage I outcome (MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,8500), breaks=c(0,2000,4000,6000, 8000))

p3.2 <-ggplot(testdf, aes(x=MM4_pred_bs1))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage I outcome (Ridge MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,8500), breaks=c(0,2000,4000,6000, 8000))

p3.3 <-ggplot(testdf, aes(x=MM5_pred_bs1))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage I outcome (Firth MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,8500), breaks=c(0,2000,4000,6000, 8000))

p3.4 <-ggplot(testdf, aes(x=MM2_pred_bs1))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage I outcome (Linear MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,8500), breaks=c(0,2000,4000,6000, 8000))

p3.5 <-ggplot(testdf, aes(x=MM7_pred_bs1))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage I outcome (RF)")+ylab("Frequency")+
  theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,8500), breaks=c(0,2000,4000,6000, 8000))

p3.6 <-ggplot(testdf, aes(x=MM8_pred_bs1))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+xlab("Estimated probability of a stage I outcome (XGBoost)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,8500), breaks=c(0,2000,4000,6000, 8000))

p3.7 <-ggplot(testdf, aes(x=MM9_pred_bs1))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage I outcome (NN)")+ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,8500), breaks=c(0,2000,4000,6000, 8000))


p3.8 <-ggplot(testdf, aes(x=MM10_pred_bs1))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+xlab("Estimated probability of a stage I outcome (SVM)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,8500), breaks=c(0,2000,4000,6000, 8000))

tiff("histograms stage 1 all models with ca125 6.tiff", width = 26, height = 31, units = "cm", res = 300)
figure <- ggarrange(p3.1, p3.2, p3.3, p3.4, p3.5,p3.6,p3.7,p3.8, ncol=2, nrow=4, common.legend=TRUE, legend="bottom", align="hv")
figure
dev.off()

#stage 2
p4.1 <-ggplot(testdf, aes(x=MM1_pred_bs2))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage II-IV outcome (MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,4000), breaks=c(0,2000,4000))

p4.2 <-ggplot(testdf, aes(x=MM4_pred_bs2))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage II-IV outcome (Ridge MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,4000), breaks=c(0,2000,4000))

p4.3 <-ggplot(testdf, aes(x=MM5_pred_bs2))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage II-IV outcome (Firth MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,4000), breaks=c(0,2000,4000))

p4.4 <-ggplot(testdf, aes(x=MM2_pred_bs2))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage II-IV outcome (Linear MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,4000), breaks=c(0,2000,4000))

p4.5 <-ggplot(testdf, aes(x=MM7_pred_bs2))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage II-IV outcome (RF)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,4000), breaks=c(0,2000,4000))

p4.6 <-ggplot(testdf, aes(x=MM8_pred_bs2))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage II-IV outcome (XGBoost)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,4000), breaks=c(0,2000,4000))

p4.7 <-ggplot(testdf, aes(x=MM9_pred_bs2))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage II-IV outcome (NN)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,4000), breaks=c(0,2000,4000))

p4.8 <-ggplot(testdf, aes(x=MM10_pred_bs2))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a stage II-IV outcome (SVM)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,4000), breaks=c(0,2000,4000))

tiff("histograms stage 2 all models with ca125 2.tiff", width = 25, height = 31, units = "cm", res = 300)
figure <- ggarrange(p4.1, p4.2, p4.3, p4.4, p4.5,p4.6,p4.7,p4.8, ncol=2, nrow=4, common.legend=TRUE, legend="bottom", align="hv")
figure
dev.off()

#metastatic
#stage 1
p5.1 <-ggplot(testdf, aes(x=MM1_pred_mt))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a secondary metastatic outcome (MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,5000), breaks=c(0,2500,5000))

p5.2 <-ggplot(testdf, aes(x=MM4_pred_mt))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a secondary metastatic outcome (Ridge MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,5000), breaks=c(0,2500,5000))

p5.3 <-ggplot(testdf, aes(x=MM5_pred_mt))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a secondary metastatic outcome (Firth MLR)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,5000), breaks=c(0,2500,5000))

p5.4 <-ggplot(testdf, aes(x=MM2_pred_mt))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a secondary metastatic outcome (Linear MLR)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,5000), breaks=c(0,2500,5000))

p5.5 <-ggplot(testdf, aes(x=MM7_pred_mt))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a secondary metastatic outcome (RF)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,5000), breaks=c(0,2500,5000))

p5.6 <-ggplot(testdf, aes(x=MM8_pred_mt))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a secondary metastatic outcome (XGBoost)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,5000), breaks=c(0,2500,5000))

p5.7 <-ggplot(testdf, aes(x=MM9_pred_mt))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a secondary metastatic outcome (NN)")+ylab("Frequency")+theme_minimal()+
  scale_x_continuous(breaks=c(0,0.5,1),labels = c("0", "0.5", "1"), limits=c(0,1))+
  scale_y_continuous(limits=c(0,5000), breaks=c(0,2500,5000))

p5.8 <-ggplot(testdf, aes(x=MM10_pred_mt))+geom_histogram(color="black", fill="white")+facet_grid(cols=vars(outcome55))+
  xlab("Estimated probability of a secondary metastatic outcome (SVM)")+
  ylab("Frequency")+theme_minimal()+scale_x_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1"), limits=c(0,1))+theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  scale_y_continuous(limits=c(0,5000), breaks=c(0,2500,5000))

tiff("histograms metastatic all models with ca125.tiff", width = 27, height = 31, units = "cm", res = 300)
figure <- ggarrange(p5.1, p5.2, p5.3, p5.4, p5.5,p5.6,p5.7,p5.8, ncol=2, nrow=4, common.legend=TRUE, legend="bottom", align="hv")
figure
dev.off()

########################
##APPARENT PERFORMANCE##
########################

setwd("~/IOTA/apparent performance plots")

#Merge the smallest centers into 1
levels(traindf$Center)[levels(traindf$Center)=="FIT"] <-"Other"
levels(traindf$Center)[levels(traindf$Center)=="GIT"] <-"Other"
levels(traindf$Center)[levels(traindf$Center)=="OCA"] <-"Other"
levels(traindf$Center)[levels(traindf$Center)=="VIT"] <-"Other"


#BENIGN VS OVERALL MALIGNANT

#AUC

#MLR

apparentmm1auc <- AUCimp.IOTA(pred=MM1_mal, outcome=outcome1, center=Center, imp=imp, data=traindf)

apparentmm1auc$Plot
apparentmm1auc$Plot[24, "RRauc"] <- "   "
apparentmm1auc$Plot[25, "RRauc"] <- "   "
apparentmm1auc$Plot[27, "RRauc"] <- "        (0.90 to 0.96)"
apparentmm1auc$Plot[27, "RRprev"] <- "   "
apparentmm1auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven", "Rome", "Malmö", "Genk", "Monza", "Prague", 
                                       "Bologna","Milan 1 ", "Lublin", "Cagliari","Milan 3", "Stockholm", "London", "Naples",
                                       "Paris","Lund", "Beijing","Other", "Maurepas","Udin","Barcelona","", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")


fpDrawBarCI <- function (lower_limit, estimate, upper_limit, size, col, y.offset = 0.5, ...){ 
  size <- ifelse(is.unit(size), convertUnit(size, unitTo = "npc", valueOnly = TRUE), size) * 0.9
  grid.polygon(x = unit(c(lower_limit, upper_limit, upper_limit, 
                          lower_limit), "native"), 
               y = c(0.4, 0.4, 0.6, 0.6), 
               gp = gpar(fill = col, col = col))
}

tiff("forest plot binary auroc apparent performance with ca125 .tiff", width = 25, height = 19, units = "cm", res = 300)

forestplot(apparentmm1auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentmm1auc$dataPlot$AUC,
           lower = apparentmm1auc$dataPlot$LL,
           upper = apparentmm1auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentmm1auc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUROC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentmm1auc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentmm1auc$Performance$AUC)
dev.off()

#LINEAR MLR
apparentmm2auc <- AUCimp.IOTA(pred=MM2_mal, outcome=outcome1, center=Center, imp=imp, data=traindf)


apparentmm2auc$Plot
apparentmm2auc$Plot[24, "RRauc"] <- "   "
apparentmm2auc$Plot[25, "RRauc"] <- "   "
apparentmm2auc$Plot[27, "RRauc"] <- "        (0.90 to 0.96)"
apparentmm2auc$Plot[27, "RRprev"] <- "   "
apparentmm2auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                       "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France","Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain","", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")




tiff("forest plot binary auroc apparent performance with ca125 linear mlr.tiff", width = 25, height = 19, units = "cm", res = 300)

forestplot(apparentmm2auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentmm2auc$dataPlot$AUC,
           lower = apparentmm2auc$dataPlot$LL,
           upper = apparentmm2auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentmm2auc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentmm2auc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentmm2auc$Performance$AUC)
dev.off()


#Ridge LR

apparentmm4auc <- AUCimp.IOTA(pred=MM4_mal, outcome=outcome1, center=Center, imp=imp, data=traindf)

apparentmm4auc$Plot
apparentmm4auc$Plot[24, "RRauc"] <- "   "
apparentmm4auc$Plot[25, "RRauc"] <- "   "
apparentmm4auc$Plot[27, "RRauc"] <- "        (0.90 to 0.96)"
apparentmm4auc$Plot[27, "RRprev"] <- "   "
apparentmm4auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                       "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France","Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain","", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")




tiff("forest plot binary auroc apparent performance with ca125 ridge mlr.tiff", width = 25, height = 19, units = "cm", res = 300)

forestplot(apparentmm4auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentmm4auc$dataPlot$AUC,
           lower = apparentmm4auc$dataPlot$LL,
           upper = apparentmm4auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentmm4auc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentmm4auc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentmm4auc$Performance$AUC)
dev.off()

#FIRTH LR

apparentmm5auc <- AUCimp.IOTA(pred=MM5_mal, outcome=outcome1, center=Center, imp=imp, data=traindf)


apparentmm5auc$Plot
apparentmm5auc$Plot[24, "RRauc"] <- "   "
apparentmm5auc$Plot[25, "RRauc"] <- "   "
apparentmm5auc$Plot[27, "RRauc"] <- "        (0.90 to 0.96)"
apparentmm5auc$Plot[27, "RRprev"] <- "   "
apparentmm5auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                       "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France","Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain","", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")



tiff("forest plot binary auroc apparent performance with ca125 firth mlr.tiff", width = 25, height = 19, units = "cm", res = 300)

forestplot(apparentmm5auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentmm5auc$dataPlot$AUC,
           lower = apparentmm5auc$dataPlot$LL,
           upper = apparentmm5auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentmm5auc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentmm5auc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentmm5auc$Performance$AUC)
dev.off()

#RANDOM FOREST

apparentmm7auc <- AUCimp.IOTA(pred=MM7_mal, outcome=outcome1, center=Center, imp=imp, data=traindf)


apparentmm7auc$Plot
apparentmm7auc$Plot[24, "RRauc"] <- "   "
apparentmm7auc$Plot[25, "RRauc"] <- "   "
apparentmm7auc$Plot[27, "RRauc"] <- "        (0.98 to 0.99)"
apparentmm7auc$Plot[27, "RRprev"] <- "   "
apparentmm7auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                       "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France","Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain","", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")



tiff("forest plot binary auroc apparent performance with ca125 rf.tiff", width = 25, height = 19, units = "cm", res = 300)

forestplot(apparentmm7auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentmm7auc$dataPlot$AUC,
           lower = apparentmm7auc$dataPlot$LL,
           upper = apparentmm7auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentmm7auc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentmm7auc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentmm7auc$Performance$AUC)
dev.off()

# XGBOOST

apparentmm8auc <- AUCimp.IOTA(pred=MM8_mal, outcome=outcome1, center=Center, imp=imp, data=traindf)


apparentmm8auc$Plot
apparentmm8auc$Plot[24, "RRauc"] <- "   "
apparentmm8auc$Plot[25, "RRauc"] <- "   "
apparentmm8auc$Plot[27, "RRauc"] <- "        (0.93 to 0.98)"
apparentmm8auc$Plot[27, "RRprev"] <- "   "
apparentmm8auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                       "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France","Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain","", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")



tiff("forest plot binary auroc apparent performance with ca125 xg.tiff", width = 25, height = 19, units = "cm", res = 300)

forestplot(apparentmm8auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentmm8auc$dataPlot$AUC,
           lower = apparentmm8auc$dataPlot$LL,
           upper = apparentmm8auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentmm8auc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentmm8auc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentmm8auc$Performance$AUC)
dev.off()

#NEURAL NETWORK

apparentmm9auc <- AUCimp.IOTA(pred=MM9_mal, outcome=outcome1, center=Center, imp=imp, data=traindf)


apparentmm9auc$Plot
apparentmm9auc$Plot[24, "RRauc"] <- "   "
apparentmm9auc$Plot[25, "RRauc"] <- "   "
apparentmm9auc$Plot[27, "RRauc"] <- "        (0.91 to 0.97)"
apparentmm9auc$Plot[27, "RRprev"] <- "   "
apparentmm9auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                       "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France","Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain","", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")



tiff("forest plot binary auroc apparent performance with ca nn.tiff", width = 25, height = 19, units = "cm", res = 300)

forestplot(apparentmm9auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentmm9auc$dataPlot$AUC,
           lower = apparentmm9auc$dataPlot$LL,
           upper = apparentmm9auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentmm9auc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentmm9auc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentmm9auc$Performance$AUC)
dev.off()

#SVM

apparentmm10auc <- AUCimp.IOTA(pred=MM10_mal, outcome=outcome1, center=Center, imp=imp, data=traindf)


apparentmm10auc$Plot
apparentmm10auc$Plot[24, "RRauc"] <- "   "
apparentmm10auc$Plot[25, "RRauc"] <- "   "
apparentmm10auc$Plot[27, "RRauc"] <- "        (0.88 to 0.96)"
apparentmm10auc$Plot[27, "RRprev"] <- "   "
apparentmm10auc$Plot[, "RRcenter"] <- c("Centre", "", "Leuven, Belgium", "Rome, Italy", "Malmö, Sweden", "Genk, Belgium", "Monza, Italy", "Prague, Czech Republic", 
                                        "Bologna, Italy","Milan 1, Italy ", "Lublin, Poland", "Cagliari, Italy","Milan 3, Italy", "Stockholm, Sweden", "London, UK", "Naples, Italy","Paris, France","Lund, Sweden", "Beijing, China","Other", "Maurepas, France","Udin, Italy","Barcelona, Spain","", "Meta-analysis", "AUC (95 % CI)", "95 % Prediction Interval")




tiff("forest plot binary auroc apparent performance with ca svm.tiff", width = 25, height = 19, units = "cm", res = 300)

forestplot(apparentmm10auc$Plot,
           align = c("l", "c", "c"),
           mean = apparentmm10auc$dataPlot$AUC,
           lower = apparentmm10auc$dataPlot$LL,
           upper = apparentmm10auc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(apparentmm10auc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(apparentmm10auc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = apparentmm10auc$Performance$AUC)
dev.off()

#Forest plot for every model

NA.forest <- apparentmm1auc$Performance[1,]
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, apparentmm1auc$Performance[1,],apparentmm4auc$Performance[1,], apparentmm5auc$Performance[1,], apparentmm2auc$Performance[1,], 
                     apparentmm7auc$Performance[1,], apparentmm8auc$Performance[1,], apparentmm9auc$Performance[1,], apparentmm10auc$Performance[1,])
Summary.AUC$Model <- c('', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'RF', 'XGBoost','NN', 'SVM')

Summary.PI <- rbind(NA.forest, apparentmm1auc$Performance[2,], apparentmm4auc$Performance[2,], apparentmm5auc$Performance[2,], apparentmm2auc$Performance[2,], 
                    apparentmm7auc$Performance[2,], apparentmm8auc$Performance[2,],apparentmm9auc$Performance[2,],apparentmm10auc$Performance[2,])
Summary.PI$Model <- c('', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'RF', 'XGBoost','NN', 'SVM')

Summary.AUC
tabletext <- cbind(
  c('Model', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'Random forest', 'XGBoost','Neural network', 'Support vector machine'),
  c('AUROC (95% CI)', 
    paste(format(round(apparentmm1auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentmm1auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentmm1auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentmm4auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentmm4auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentmm4auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentmm5auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentmm5auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentmm5auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentmm2auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentmm2auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentmm2auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentmm7auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentmm7auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentmm7auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentmm8auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentmm8auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentmm8auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentmm9auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentmm9auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentmm9auc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(apparentmm10auc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(apparentmm10auc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(apparentmm10auc$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(apparentmm1auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentmm1auc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(apparentmm4auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentmm4auc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(apparentmm5auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentmm5auc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(apparentmm2auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentmm2auc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(apparentmm7auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentmm7auc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(apparentmm8auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentmm8auc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(apparentmm9auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentmm9auc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(apparentmm10auc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(apparentmm10auc$Performance$UL[2], 2), nsmall = 2), ")"))
)


tiff("summary plot of binary auroc apparent with ca .tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE),
           xlab = "AUROC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.9, 0.95, 1), xlog = TRUE, clip = c(0.9, 1))
dev.off()



######
#PDI#
#####


#MLR


ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  PDIM1 = ests(y=tmp$outcome5, d=tmp[, c(47:51)], acc="pdi", level=0.95, method="prob", k=5)
  
  
  return(PDIM1)
}, simplify = F, USE.NAMES = T)

ResultsPDIM1APP = do.call("sapply", args = ArgzCalc)

MM1PDIapp = do.call("rbind", lapply(ResultsPDIM1APP, "[[", "value")) 
MM1PDISEapp = do.call("rbind", lapply(ResultsPDIM1APP, "[[", "se")) 


#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM1PDIapp)
WithinVar <- mean((MM1PDISEapp)^2)
BetweenVar <- var(MM1PDIapp)
PooledVar <- WithinVarMM1 + BetweenVarMM1 + BetweenVarMM1/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined

#LINEAR MLR
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  PDIM2 = ests(y=tmp$outcome5, d=tmp[, c(53:57)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM2)
}, simplify = F, USE.NAMES = T)

ResultsPDIM2APP = do.call("sapply", args = ArgzCalc)

MM2PDIapp = do.call("rbind", lapply(ResultsPDIM2APP, "[[", "value")) 
MM2PDISEapp = do.call("rbind", lapply(ResultsPDIM2APP, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM2PDIapp)
WithinVar <- mean((MM2PDISEapp)^2)
BetweenVar <- var(MM2PDIapp)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined


#RIDGE LR
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  PDIM4 = ests(y=tmp$outcome5, d=tmp[, c(65:69)], acc="pdi", level=0.95, method="prob", k=5) #change here
  
  
  return(PDIM4)
}, simplify = F, USE.NAMES = T)

ResultsPDIM4APP = do.call("sapply", args = ArgzCalc)

MM4PDIapp = do.call("rbind", lapply(ResultsPDIM4APP, "[[", "value")) 
MM4PDISEapp = do.call("rbind", lapply(ResultsPDIM4APP, "[[", "se")) 


#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM4PDIapp)
WithinVar <- mean((MM4PDISEapp)^2)
BetweenVar <- var(MM4PDIapp)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined

#FIRTH LR
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  PDIM5 = ests(y=tmp$outcome5, d=tmp[, c(102:106)], acc="pdi", level=0.95, method="prob", k=5) #change here
  
  
  return(PDIM5)
}, simplify = F, USE.NAMES = T)

ResultsPDIM5APP = do.call("sapply", args = ArgzCalc)

MM5PDIapp = do.call("rbind", lapply(ResultsPDIM5APP, "[[", "value")) 
MM5PDISEapp = do.call("rbind", lapply(ResultsPDIM5APP, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM5PDIapp)
WithinVar <- mean((MM5PDISEapp)^2)
BetweenVar <- var(MM5PDIapp)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined


#RF
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  PDIM7 = ests(y=tmp$outcome5, d=tmp[, c(48:52)], acc="pdi", level=0.95, method="prob", k=5)
  
  
  return(PDIM7)
}, simplify = F, USE.NAMES = T)

ResultsPDIM7APP = do.call("sapply", args = ArgzCalc)

MM7PDIapp = do.call("rbind", lapply(ResultsPDIM7APP, "[[", "value")) 
MM7PDISEapp = do.call("rbind", lapply(ResultsPDIM7APP, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM7PDIapp)
WithinVar <- mean((MM7PDISEapp)^2)
BetweenVar <- var(MM7PDIapp)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined


#XGBOOST
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  PDIM8 = ests(y=tmp$outcome5, d=tmp[, c(96:100)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM8)
}, simplify = F, USE.NAMES = T)

ResultsPDIM8APP = do.call("sapply", args = ArgzCalc)

MM8PDIapp = do.call("rbind", lapply(ResultsPDIM8APP, "[[", "value")) 
MM8PDISEapp = do.call("rbind", lapply(ResultsPDIM8APP, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM8PDIapp)
WithinVar <- mean((MM8PDISEapp)^2)
BetweenVar <- var(MM8PDIapp)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined


#NN

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  PDIM9 = ests(y=tmp$outcome5, d=tmp[, c(77:81)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM9)
}, simplify = F, USE.NAMES = T)

ResultsPDIM9APP = do.call("sapply", args = ArgzCalc)

MM9PDIapp = do.call("rbind", lapply(ResultsPDIM9APP, "[[", "value")) 
MM9PDISEapp = do.call("rbind", lapply(ResultsPDIM9APP, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM9PDIapp)
WithinVar <- mean((MM9PDISEapp)^2)
BetweenVar <- var(MM9PDIapp)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined

#SVM
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  PDIM10 = ests(y=tmp$outcome5, d=tmp[, c(83:87)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM10)
}, simplify = F, USE.NAMES = T)

ResultsPDIM10APP = do.call("sapply", args = ArgzCalc)

MM10PDIapp = do.call("rbind", lapply(ResultsPDIM10APP, "[[", "value")) 
MM10PDISEapp = do.call("rbind", lapply(ResultsPDIM10APP, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM10PDIapp)
WithinVar <- mean((MM10PDISEapp)^2)
BetweenVar <- var(MM10PDIapp)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined


#c statistics for every pair of outcomes

#MLR

#pair 1 benign vs borderline
traindf$bvsb <- traindf$MM1_pred_be / (traindf$MM1_pred_be + traindf$MM1_pred_bo)
AUC.imp(pred = traindf$bvsb, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
traindf$bvss1 <- traindf$MM1_pred_be / (traindf$MM1_pred_be +traindf$MM1_pred_bs1)
AUC.imp(pred = traindf$bvss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
traindf$bvss2 <- traindf$MM1_pred_be / (traindf$MM1_pred_be +traindf$MM1_pred_bs2)
AUC.imp(pred = traindf$bvss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
traindf$bvsmt <- traindf$MM1_pred_be / (traindf$MM1_pred_be +traindf$MM1_pred_mt)
AUC.imp(pred = traindf$bvsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
traindf$bovss1 <- traindf$MM1_pred_bo / (traindf$MM1_pred_bo +traindf$MM1_pred_bs1)
AUC.imp(pred = traindf$bovss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
traindf$bovss2 <- traindf$MM1_pred_bo / (traindf$MM1_pred_bo +traindf$MM1_pred_bs2)
AUC.imp(pred = traindf$bovss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
traindf$bovsmt <- traindf$MM1_pred_bo / (traindf$MM1_pred_bo +traindf$MM1_pred_mt)
AUC.imp(pred = traindf$bovsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
traindf$s1vss2 <- traindf$MM1_pred_bs1 / (traindf$MM1_pred_bs1 +traindf$MM1_pred_bs2)
AUC.imp(pred = traindf$s1vss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
traindf$s1vsmt <- traindf$MM1_pred_bs1 / (traindf$MM1_pred_bs1 +traindf$MM1_pred_mt)
AUC.imp(pred = traindf$s1vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
traindf$s2vsmt <- traindf$MM1_pred_bs2 / (traindf$MM1_pred_bs2 +traindf$MM1_pred_mt)
AUC.imp(pred = traindf$s2vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 4 | traindf$outcome5 == 5,])$Performance

#LINEAR MLR

#pair 1 benign vs borderline
traindf$bvsbm2 <- traindf$MM2_pred_be / (traindf$MM2_pred_be + traindf$MM2_pred_bo)
AUC.imp(pred = traindf$bvsbm2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
traindf$bvss1m2 <- traindf$MM2_pred_be / (traindf$MM2_pred_be +traindf$MM2_pred_bs1)
AUC.imp(pred = traindf$bvss1m2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
traindf$bvss2m2 <- traindf$MM2_pred_be / (traindf$MM2_pred_be +traindf$MM2_pred_bs2)
AUC.imp(pred = traindf$bvss2m2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
traindf$bvsmtm2 <- traindf$MM2_pred_be / (traindf$MM2_pred_be +traindf$MM2_pred_mt)
AUC.imp(pred = traindf$bvsmtm2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
traindf$bovss1m2 <- traindf$MM2_pred_bo / (traindf$MM2_pred_bo +traindf$MM2_pred_bs1)
AUC.imp(pred = traindf$bovss1m2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
traindf$bovss2m2 <- traindf$MM2_pred_bo / (traindf$MM2_pred_bo +traindf$MM2_pred_bs2)
AUC.imp(pred = traindf$bovss2m2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
traindf$bovsmtm2 <- traindf$MM2_pred_bo / (traindf$MM2_pred_bo +traindf$MM2_pred_mt)
AUC.imp(pred = traindf$bovsmtm2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
traindf$s1vss2m2 <- traindf$MM2_pred_bs1 / (traindf$MM2_pred_bs1 +traindf$MM2_pred_bs2)
AUC.imp(pred = traindf$s1vss2m2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
traindf$s1vsmtm2 <- traindf$MM2_pred_bs1 / (traindf$MM2_pred_bs1 +traindf$MM2_pred_mt)
AUC.imp(pred = traindf$s1vsmtm2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
traindf$s2vsmtm2 <- traindf$MM2_pred_bs2 / (traindf$MM2_pred_bs2 +traindf$MM2_pred_mt)
AUC.imp(pred = traindf$s2vsmtm2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 4 | traindf$outcome5 == 5,])$Performance


#RIDGE
#pair 1 benign vs borderline
traindf$bvsb <- traindf$MM4_pred_be / (traindf$MM4_pred_be + traindf$MM4_pred_bo)
AUC.imp(pred = traindf$bvsb, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
traindf$bvss1 <- traindf$MM4_pred_be / (traindf$MM4_pred_be +traindf$MM4_pred_bs1)
AUC.imp(pred = traindf$bvss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
traindf$bvss2 <- traindf$MM4_pred_be / (traindf$MM4_pred_be +traindf$MM4_pred_bs2)
AUC.imp(pred = traindf$bvss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
traindf$bvsmt <- traindf$MM4_pred_be / (traindf$MM4_pred_be +traindf$MM4_pred_mt)
AUC.imp(pred = traindf$bvsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
traindf$bovss1 <- traindf$MM4_pred_bo / (traindf$MM4_pred_bo +traindf$MM4_pred_bs1)
AUC.imp(pred = traindf$bovss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
traindf$bovss2 <- traindf$MM4_pred_bo / (traindf$MM4_pred_bo +traindf$MM4_pred_bs2)
AUC.imp(pred = traindf$bovss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
traindf$bovsmt <- traindf$MM4_pred_bo / (traindf$MM4_pred_bo +traindf$MM4_pred_mt)
AUC.imp(pred = traindf$bovsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
traindf$s1vss2 <- traindf$MM4_pred_bs1 / (traindf$MM4_pred_bs1 +traindf$MM4_pred_bs2)
AUC.imp(pred = traindf$s1vss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
traindf$s1vsmt <- traindf$MM4_pred_bs1 / (traindf$MM4_pred_bs1 +traindf$MM4_pred_mt)
AUC.imp(pred = traindf$s1vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
traindf$s2vsmt <- traindf$MM4_pred_bs2 / (traindf$MM4_pred_bs2 +traindf$MM4_pred_mt)
AUC.imp(pred = traindf$s2vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 4 | traindf$outcome5 == 5,])$Performance

#FIRTH
#pair 1 benign vs borderline
traindf$bvsb <- traindf$MM5_pred_be / (traindf$MM5_pred_be + traindf$MM5_pred_bo)
AUC.imp(pred = traindf$bvsb, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
traindf$bvss1 <- traindf$MM5_pred_be / (traindf$MM5_pred_be +traindf$MM5_pred_bs1)
AUC.imp(pred = traindf$bvss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
traindf$bvss2 <- traindf$MM5_pred_be / (traindf$MM5_pred_be +traindf$MM5_pred_bs2)
AUC.imp(pred = traindf$bvss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
traindf$bvsmt <- traindf$MM5_pred_be / (traindf$MM5_pred_be +traindf$MM5_pred_mt)
AUC.imp(pred = traindf$bvsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
traindf$bovss1 <- traindf$MM5_pred_bo / (traindf$MM5_pred_bo +traindf$MM5_pred_bs1)
AUC.imp(pred = traindf$bovss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
traindf$bovss2 <- traindf$MM5_pred_bo / (traindf$MM5_pred_bo +traindf$MM5_pred_bs2)
AUC.imp(pred = traindf$bovss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
traindf$bovsmt <- traindf$MM5_pred_bo / (traindf$MM5_pred_bo +traindf$MM5_pred_mt)
AUC.imp(pred = traindf$bovsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
traindf$s1vss2 <- traindf$MM5_pred_bs1 / (traindf$MM5_pred_bs1 +traindf$MM5_pred_bs2)
AUC.imp(pred = traindf$s1vss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 4,])$Performance


#pair 9: stage I vs secondary metastatic
traindf$s1vsmt <- traindf$MM5_pred_bs1 / (traindf$MM5_pred_bs1 +traindf$MM5_pred_mt)
AUC.imp(pred = traindf$s1vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 5,])$Performance


#pair 10: stage II-IV vs secondary metastatic
traindf$s2vsmt <- traindf$MM5_pred_bs2 / (traindf$MM5_pred_bs2 +traindf$MM5_pred_mt)
AUC.imp(pred = traindf$s2vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 4 | traindf$outcome5 == 5,])$Performance


#rf
#pair 1 benign vs borderline
traindf$bvsb <- traindf$MM7_pred_be / (traindf$MM7_pred_be + traindf$MM7_pred_bo)
AUC.imp(pred = traindf$bvsb, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
traindf$bvss1 <- traindf$MM7_pred_be / (traindf$MM7_pred_be +traindf$MM7_pred_bs1)
AUC.imp(pred = traindf$bvss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
traindf$bvss2 <- traindf$MM7_pred_be / (traindf$MM7_pred_be +traindf$MM7_pred_bs2)
AUC.imp(pred = traindf$bvss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
traindf$bvsmt <- traindf$MM7_pred_be / (traindf$MM7_pred_be +traindf$MM7_pred_mt)
AUC.imp(pred = traindf$bvss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
traindf$bovss1 <- traindf$MM7_pred_bo / (traindf$MM7_pred_bo +traindf$MM7_pred_bs1)
AUC.imp(pred = traindf$bovss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
traindf$bovss2 <- traindf$MM7_pred_bo / (traindf$MM7_pred_bo +traindf$MM7_pred_bs2)
AUC.imp(pred = traindf$bovss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
traindf$bovsmt <- traindf$MM7_pred_bo / (traindf$MM7_pred_bo +traindf$MM7_pred_mt)
AUC.imp(pred = traindf$bovsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
traindf$s1vss2 <- traindf$MM7_pred_bs1 / (traindf$MM7_pred_bs1 +traindf$MM7_pred_bs2)
AUC.imp(pred = traindf$s1vss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
traindf$s1vsmt <- traindf$MM7_pred_bs1 / (traindf$MM7_pred_bs1 +traindf$MM7_pred_mt)
AUC.imp(pred = traindf$s1vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
traindf$s2vsmt <- traindf$MM7_pred_bs2 / (traindf$MM7_pred_bs2 +traindf$MM7_pred_mt)
AUC.imp(pred = traindf$s2vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 4 | traindf$outcome5 == 5,])$Performance

#xgboost
#pair 1 benign vs borderline
traindf$bvsb <- traindf$MM8_pred_be / (traindf$MM8_pred_be + traindf$MM8_pred_bo)
AUC.imp(pred = traindf$bvsb, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
traindf$bvss1 <- traindf$MM8_pred_be / (traindf$MM8_pred_be +traindf$MM8_pred_bs1)
AUC.imp(pred = traindf$bvss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
traindf$bvss2 <- traindf$MM8_pred_be / (traindf$MM8_pred_be +traindf$MM8_pred_bs2)
AUC.imp(pred = traindf$bvss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
traindf$bvsmt <- traindf$MM8_pred_be / (traindf$MM8_pred_be +traindf$MM8_pred_mt)
AUC.imp(pred = traindf$bvsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
traindf$bovss1 <- traindf$MM8_pred_bo / (traindf$MM8_pred_bo +traindf$MM8_pred_bs1)
AUC.imp(pred = traindf$bovss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
traindf$bovss2 <- traindf$MM8_pred_bo / (traindf$MM8_pred_bo +traindf$MM8_pred_bs2)
AUC.imp(pred = traindf$bovss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
traindf$bovsmt <- traindf$MM8_pred_bo / (traindf$MM8_pred_bo +traindf$MM8_pred_mt)
AUC.imp(pred = traindf$bovsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
traindf$s1vss2 <- traindf$MM8_pred_bs1 / (traindf$MM8_pred_bs1 +traindf$MM8_pred_bs2)
AUC.imp(pred = traindf$s1vss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
traindf$s1vsmt <- traindf$MM8_pred_bs1 / (traindf$MM8_pred_bs1 +traindf$MM8_pred_mt)
AUC.imp(pred = traindf$s1vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
traindf$s2vsmt <- traindf$MM8_pred_bs2 / (traindf$MM8_pred_bs2 +traindf$MM8_pred_mt)
AUC.imp(pred = traindf$s2vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 4 | traindf$outcome5 == 5,])$Performance

#nn
#pair 1 benign vs borderline
traindf$bvsb <- traindf$MM9_pred_be / (traindf$MM9_pred_be + traindf$MM9_pred_bo)
AUC.imp(pred = traindf$bvsb, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
traindf$bvss1 <- traindf$MM9_pred_be / (traindf$MM9_pred_be +traindf$MM9_pred_bs1)
AUC.imp(pred = traindf$bvss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
traindf$bvss2 <- traindf$MM9_pred_be / (traindf$MM9_pred_be +traindf$MM9_pred_bs2)
AUC.imp(pred = traindf$bvss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
traindf$bvsmt <- traindf$MM9_pred_be / (traindf$MM9_pred_be +traindf$MM9_pred_mt)
AUC.imp(pred = traindf$bvsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
traindf$bovss1 <- traindf$MM9_pred_bo / (traindf$MM9_pred_bo +traindf$MM9_pred_bs1)
AUC.imp(pred = traindf$bovss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
traindf$bovss2 <- traindf$MM9_pred_bo / (traindf$MM9_pred_bo +traindf$MM9_pred_bs2)
AUC.imp(pred = traindf$bovss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
traindf$bovsmt <- traindf$MM9_pred_bo / (traindf$MM9_pred_bo +traindf$MM9_pred_mt)
AUC.imp(pred = traindf$bovsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
traindf$s1vss2 <- traindf$MM9_pred_bs1 / (traindf$MM9_pred_bs1 +traindf$MM9_pred_bs2)
AUC.imp(pred = traindf$s1vss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
traindf$s1vsmt <- traindf$MM9_pred_bs1 / (traindf$MM9_pred_bs1 +traindf$MM9_pred_mt)
AUC.imp(pred = traindf$s1vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
traindf$s2vsmt <- traindf$MM9_pred_bs2 / (traindf$MM9_pred_bs2 +traindf$MM9_pred_mt)
AUC.imp(pred = traindf$s2vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 4 | traindf$outcome5 == 5,])$Performance

#svm
#pair 1 benign vs borderline
traindf$bvsb <- traindf$MM10_pred_be / (traindf$MM10_pred_be + traindf$MM10_pred_bo)
AUC.imp(pred = traindf$bvsb, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
traindf$bvss1 <- traindf$MM10_pred_be / (traindf$MM10_pred_be +traindf$MM10_pred_bs1)
AUC.imp(pred = traindf$bvss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
traindf$bvss2 <- traindf$MM10_pred_be / (traindf$MM10_pred_be +traindf$MM10_pred_bs2)
AUC.imp(pred = traindf$bvss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
traindf$bvsmt <- traindf$MM10_pred_be / (traindf$MM10_pred_be +traindf$MM10_pred_mt)
AUC.imp(pred = traindf$bvsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 1 | traindf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
traindf$bovss1 <- traindf$MM10_pred_bo / (traindf$MM10_pred_bo +traindf$MM10_pred_bs1)
AUC.imp(pred = traindf$bovss1, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
traindf$bovss2 <- traindf$MM10_pred_bo / (traindf$MM10_pred_bo +traindf$MM10_pred_bs2)
AUC.imp(pred = traindf$bovss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
traindf$bovsmt <- traindf$MM10_pred_bo / (traindf$MM10_pred_bo +traindf$MM10_pred_mt)
AUC.imp(pred = traindf$bovsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 2 | traindf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
traindf$s1vss2 <- traindf$MM10_pred_bs1 / (traindf$MM10_pred_bs1 +traindf$MM10_pred_bs2)
AUC.imp(pred = traindf$s1vss2, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
traindf$s1vsmt <- traindf$MM10_pred_bs1 / (traindf$MM10_pred_bs1 +traindf$MM10_pred_mt)
AUC.imp(pred = traindf$s1vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 3 | traindf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
traindf$s2vsmt <- traindf$MM10_pred_bs2 / (traindf$MM10_pred_bs2 +traindf$MM10_pred_mt)
AUC.imp(pred = traindf$s2vsmt, outcome = traindf$outcome5, imp = traindf$imp, data = traindf[traindf$outcome5 == 4 | traindf$outcome5 == 5,])$Performance


###########################
####EXTERNAL VALIDATION####
###########################

setwd("~/Travail/External validation ML update")

levels(testdf$Center)[levels(testdf$Center)=="IUK"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="MIT"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="FLI"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="NUK"] <-"Other"


#AUC FOR BENIGN VS MALIGNANT

#MLR
mm1_extauc <- AUCimp.IOTA(pred=MM1_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)


mm1_extauc$Plot
mm1_extauc$Plot[17, "RRauc"] <- "   "
mm1_extauc$Plot[18, "RRauc"] <- "   "
mm1_extauc$Plot[20, "RRauc"] <- "        (0.86 to 0.96)"
mm1_extauc$Plot[20, "RRprev"] <- "   "

mm1_extauc$Plot[, "RRcenter"] <- c("Centre", "", "Athens", "Rome", "Milan 1", "Malmö", "Stockholm", "Genk", "Leuven", 
                                   "Monza","Cagliari", "Other", "Milan 2", "Pamplona", "Trieste", "Katowice","", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")



tiff("forest plot binary auroc ext val with ca mlr.tiff", width = 25, height = 17.75, units = "cm", res = 300)

forestplot(mm1_extauc$Plot,
           align = c("l", "c", "c"),
           mean = mm1_extauc$dataPlot$AUC,
           lower = mm1_extauc$dataPlot$LL,
           upper = mm1_extauc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm1_extauc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUROC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm1_extauc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm1_extauc$Performance$AUC)

dev.off()

#LINEAR MLR

mm2_extauc <- AUCimp.IOTA(pred=MM2_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)


mm2_extauc$Plot
mm2_extauc$Plot[17, "RRauc"] <- "   "
mm2_extauc$Plot[18, "RRauc"] <- "   "
mm2_extauc$Plot[20, "RRauc"] <- "        (0.84 to 0.94)"
mm2_extauc$Plot[20, "RRprev"] <- "   "

mm2_extauc$Plot[, "RRcenter"] <- c("Centre", "", "Athens", "Rome", "Milan 1", "Malmö", "Stockholm", "Genk", "Leuven", 
                                   "Monza","Cagliari", "Other", "Milan 2", "Pamplona", "Trieste", "Katowice","", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")



tiff("forest plot binary auroc ext val with ca linear mlr.tiff", width = 25, height = 17.75, units = "cm", res = 300)

forestplot(mm2_extauc$Plot,
           align = c("l", "c", "c"),
           mean = mm2_extauc$dataPlot$AUC,
           lower = mm2_extauc$dataPlot$LL,
           upper = mm2_extauc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm2_extauc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUROC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm2_extauc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm2_extauc$Performance$AUC)

dev.off()


#RIDGE
mm4_extauc <- AUCimp.IOTA(pred=MM4_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)

mm4_extauc$Plot
mm4_extauc$Plot[17, "RRauc"] <- "   "
mm4_extauc$Plot[18, "RRauc"] <- "   "
mm4_extauc$Plot[20, "RRauc"] <- "        (0.86 to 0.95)"
mm4_extauc$Plot[20, "RRprev"] <- "   "

mm4_extauc$Plot[, "RRcenter"] <- c("Centre", "", "Athens", "Rome", "Milan 1", "Malmö", "Stockholm", "Genk", "Leuven", 
                                   "Monza","Cagliari", "Other", "Milan 2", "Pamplona", "Trieste", "Katowice","", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")



tiff("forest plot binary auroc ext val with ca ridge mlr 2.tiff", width = 25, height = 17.75, units = "cm", res = 300)

forestplot(mm4_extauc$Plot,
           align = c("l", "c", "c"),
           mean = mm4_extauc$dataPlot$AUC,
           lower = mm4_extauc$dataPlot$LL,
           upper = mm4_extauc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm4_extauc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUROC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm4_extauc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm4_extauc$Performance$AUC)

dev.off()

#FIRTH
mm5_extauc <- AUCimp.IOTA(pred=MM5_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)


mm5_extauc$Plot
mm5_extauc$Plot[17, "RRauc"] <- "   "
mm5_extauc$Plot[18, "RRauc"] <- "   "
mm5_extauc$Plot[20, "RRauc"] <- "        (0.86 to 0.96)"
mm5_extauc$Plot[20, "RRprev"] <- "   "

mm5_extauc$Plot[, "RRcenter"] <- c("Centre", "", "Athens", "Rome", "Milan 1", "Malmö", "Stockholm", "Genk", "Leuven", 
                                   "Monza","Cagliari", "Other", "Milan 2", "Pamplona", "Trieste", "Katowice","", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")



tiff("forest plot binary auroc ext val with ca firth mlr.tiff", width = 25, height = 17.75, units = "cm", res = 300)

forestplot(mm5_extauc$Plot,
           align = c("l", "c", "c"),
           mean = mm5_extauc$dataPlot$AUC,
           lower = mm5_extauc$dataPlot$LL,
           upper = mm5_extauc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm5_extauc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUROC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm5_extauc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm5_extauc$Performance$AUC)

dev.off()



#RF
mm7_extauc <- AUCimp.IOTA(pred=MM7_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)


mm7_extauc$Plot
mm7_extauc$Plot[17, "RRauc"] <- "   "
mm7_extauc$Plot[18, "RRauc"] <- "   "
mm7_extauc$Plot[20, "RRauc"] <- "        (0.86 to 0.96)"
mm7_extauc$Plot[20, "RRprev"] <- "   "

mm7_extauc$Plot[, "RRcenter"] <- c("Centre", "", "Athens", "Rome", "Milan 1", "Malmö", "Stockholm", "Genk", "Leuven", 
                                   "Monza","Cagliari", "Other", "Milan 2", "Pamplona", "Trieste", "Katowice","", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")



tiff("forest plot binary auroc ext val with ca rf 2.tiff", width = 25, height = 17.75, units = "cm", res = 300)

forestplot(mm7_extauc$Plot,
           align = c("l", "c", "c"),
           mean = mm7_extauc$dataPlot$AUC,
           lower = mm7_extauc$dataPlot$LL,
           upper = mm7_extauc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm7_extauc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUROC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm7_extauc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm7_extauc$Performance$AUC)

dev.off()


#XG
mm8_extauc <- AUCimp.IOTA(pred=MM8_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)

mm8_extauc$Plot
mm8_extauc$Plot[17, "RRauc"] <- "   "
mm8_extauc$Plot[18, "RRauc"] <- "   "
mm8_extauc$Plot[20, "RRauc"] <- "        (0.86 to 0.96)"
mm8_extauc$Plot[20, "RRprev"] <- "   "

mm8_extauc$Plot[, "RRcenter"] <- c("Centre", "", "Athens", "Rome", "Milan 1", "Malmö", "Stockholm", "Genk", "Leuven", 
                                   "Monza","Cagliari", "Other", "Milan 2", "Pamplona", "Trieste", "Katowice","", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")

tiff("forest plot binary auroc ext val with ca xgboost.tiff", width = 25, height = 17.75, units = "cm", res = 300)

forestplot(mm8_extauc$Plot,
           align = c("l", "c", "c"),
           mean = mm8_extauc$dataPlot$AUC,
           lower = mm8_extauc$dataPlot$LL,
           upper = mm8_extauc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm8_extauc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUROC",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm8_extauc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm8_extauc$Performance$AUC)

dev.off()

#NN
mm9_extauc <- AUCimp.IOTA(pred=MM9_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)


mm9_extauc$Plot
mm9_extauc$Plot[17, "RRauc"] <- "   " 
mm9_extauc$Plot[18, "RRauc"] <- "   "
mm9_extauc$Plot[20, "RRauc"] <- "        (0.86 to 0.96)"
mm9_extauc$Plot[20, "RRprev"] <- "   "

mm9_extauc$Plot[, "RRcenter"] <- c("Centre", "", "Athens", "Rome", "Milan 1", "Malmö", "Stockholm", "Genk", "Leuven", 
                                   "Monza","Cagliari", "Other", "Milan 2", "Pamplona", "Trieste", "Katowice","", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")

tiff("forest plot binary auroc ext val wit ca nn.tiff", width = 25, height = 17.75, units = "cm", res = 300)

forestplot(mm9_extauc$Plot,
           align = c("l", "c", "c"),
           mean = mm9_extauc$dataPlot$AUC,
           lower = mm9_extauc$dataPlot$LL,
           upper = mm9_extauc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm9_extauc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUROC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm9_extauc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm9_extauc$Performance$AUC)
dev.off()

#SVM
mm10_extauc <- AUCimp.IOTA(pred=MM10_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)


mm10_extauc$Plot
mm10_extauc$Plot[17, "RRauc"] <- "   "
mm10_extauc$Plot[18, "RRauc"] <- "   "
mm10_extauc$Plot[20, "RRauc"] <- "        (0.87 to 0.91)"
mm10_extauc$Plot[20, "RRprev"] <- "   "

mm10_extauc$Plot[, "RRcenter"] <- c("Centre", "", "Athens", "Rome", "Milan 1", "Malmö", "Stockholm", "Genk", "Leuven", 
                                    "Monza","Cagliari", "Other", "Milan 2", "Pamplona", "Trieste", "Katowice","", "Meta-analysis", "AUROC (95 % CI)", "95 % Prediction Interval")


tiff("forest plot binary auroc ext val with ca svm.tiff", width = 25, height = 17.75, units = "cm", res = 300)

forestplot(mm10_extauc$Plot,
           align = c("l", "c", "c"),
           mean = mm10_extauc$dataPlot$AUC,
           lower = mm10_extauc$dataPlot$LL,
           upper = mm10_extauc$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(mm10_extauc$IncludedCenters)), TRUE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUROC", 
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(mm10_extauc$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = mm10_extauc$Performance$AUC)

dev.off()

#Forest plot for every model

NA.forest <- mm1_extauc$Performance[1,]
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, mm1_extauc$Performance[1,],mm4_extauc$Performance[1,], mm5_extauc$Performance[1,], mm2_extauc$Performance[1,],
                     mm7_extauc$Performance[1,], mm8_extauc$Performance[1,], mm9_extauc$Performance[1,], mm10_extauc$Performance[1,])
Summary.AUC$Model <- c('', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'RF', 'XGBoost','NN', 'SVM')

Summary.PI <- rbind(NA.forest, mm1_extauc$Performance[2,], mm4_extauc$Performance[2,], mm5_extauc$Performance[2,], mm2_extauc$Performance[2,],
                    mm7_extauc$Performance[2,], mm8_extauc$Performance[2,],mm9_extauc$Performance[2,],mm10_extauc$Performance[2,])
Summary.PI$Model <- c('', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'RF', 'XGBoost','NN', 'SVM')

Summary.AUC
tabletext <- cbind(
  c('Model', 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'Random forest', 'XGBoost','Neural network', 'Support vector machine'),
  c('AUROC (95% CI)', 
    paste(format(round(mm1_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm1_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm1_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm4_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm4_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm4_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm5_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm5_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm5_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm2_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm2_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm2_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm7_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm7_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm7_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm8_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm8_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm8_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm9_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm9_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm9_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm10_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm10_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm10_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(mm1_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm1_extauc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm4_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm4_extauc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm5_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm5_extauc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm2_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm2_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm7_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm7_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm8_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm8_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm9_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm9_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm10_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm10_extauc$Performance$UL[2], 2), nsmall = 2), ")"))
)


tiff("summary plot of binary auroc ext val with ca.tiff", width = 31, height = 20, units = "cm", res = 300)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE, TRUE,  TRUE, TRUE),
           xlab = "AUROC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.85, 1))
dev.off()
#############
#PDI

#MLR

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = testdf[testdf$imp == i, ]
  PDIM1 = ests(y=tmp$outcome5, d=tmp[, c(48:52)], acc="pdi", level=0.95, method="prob", k=5)
  
  
  return(PDIM1)
}, simplify = F, USE.NAMES = T)

ResultsPDIM1 = do.call("sapply", args = ArgzCalc)

MM1PDI = do.call("rbind", lapply(ResultsPDIM1, "[[", "value")) 
MM1PDISE = do.call("rbind", lapply(ResultsPDIM1, "[[", "se")) 


#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM1PDI)
WithinVar <- mean((MM1PDISE)^2)
BetweenVar <- var(MM1PDI)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined

#LINEAR MLR
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = testdf[testdf$imp == i, ]
  PDIM2 = ests(y=tmp$outcome5, d=tmp[, c(54:58)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM2)
}, simplify = F, USE.NAMES = T)

ResultsPDIM2 = do.call("sapply", args = ArgzCalc)

MM2PDI = do.call("rbind", lapply(ResultsPDIM2, "[[", "value")) 
MM2PDISE = do.call("rbind", lapply(ResultsPDIM2, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM2PDI)
WithinVar <- mean((MM2PDISE)^2)
BetweenVar <- var(MM2PDI)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined



#RIDGE LR
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = testdf[testdf$imp == i, ]
  PDIM4 = ests(y=tmp$outcome5, d=tmp[, c(66:70)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM4)
}, simplify = F, USE.NAMES = T)

ResultsPDIM4 = do.call("sapply", args = ArgzCalc)

MM4PDI = do.call("rbind", lapply(ResultsPDIM4, "[[", "value")) 
MM4PDISE = do.call("rbind", lapply(ResultsPDIM4, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM4PDI)
WithinVar <- mean((MM4PDISE)^2)
BetweenVar <- var(MM4PDI)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined

#FIRTH LR
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = testdf[testdf$imp == i, ]
  PDIM5 = ests(y=tmp$outcome5, d=tmp[, c(102:106)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM5)
}, simplify = F, USE.NAMES = T)

ResultsPDIM5 = do.call("sapply", args = ArgzCalc)

MM5PDI = do.call("rbind", lapply(ResultsPDIM5, "[[", "value")) 
MM5PDISE = do.call("rbind", lapply(ResultsPDIM5, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM5PDI)
WithinVar <- mean((MM5PDISE)^2)
BetweenVar <- var(MM5PDI)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined


#RF
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = testdf[testdf$imp == i, ]
  PDIM7 = ests(y=tmp$outcome5, d=tmp[, c(48:52)], acc="pdi", level=0.95, method="prob", k=5) 
  
  return(PDIM7)
}, simplify = F, USE.NAMES = T)

ResultsPDIM7 = do.call("sapply", args = ArgzCalc)

MM7PDI = do.call("rbind", lapply(ResultsPDIM7, "[[", "value")) 
MM7PDISE = do.call("rbind", lapply(ResultsPDIM7, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM7PDI)
WithinVar <- mean((MM7PDISE)^2)
BetweenVar <- var(MM7PDI)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined


#XGBOOST
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = testdf[testdf$imp == i, ]
  PDIM8 = ests(y=tmp$outcome5, d=tmp[, c(84:88)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM8)
}, simplify = F, USE.NAMES = T)

ResultsPDIM8 = do.call("sapply", args = ArgzCalc)

MM8PDI = do.call("rbind", lapply(ResultsPDIM8, "[[", "value")) 
MM8PDISE = do.call("rbind", lapply(ResultsPDIM8, "[[", "se")) 


#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM8PDI)
WithinVar <- mean((MM8PDISE)^2)
BetweenVar <- var(MM8PDI)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined


#NN

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = testdf[testdf$imp == i, ]
  PDIM9 = ests(y=tmp$outcome5, d=tmp[, c(90:94)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM9)
}, simplify = F, USE.NAMES = T)

ResultsPDIM9 = do.call("sapply", args = ArgzCalc)

MM9PDI = do.call("rbind", lapply(ResultsPDIM9, "[[", "value")) 
MM9PDISE = do.call("rbind", lapply(ResultsPDIM9, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM9PDI)
WithinVar <- mean((MM9PDISE)^2)
BetweenVar <- var(MM9PDI)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined

#SVM
ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = testdf[testdf$imp == i, ]
  PDIM10 = ests(y=tmp$outcome5, d=tmp[, c(96:100)], acc="pdi", level=0.95, method="prob", k=5) 
  
  
  return(PDIM10)
}, simplify = F, USE.NAMES = T)

ResultsPDIM10 = do.call("sapply", args = ArgzCalc)

MM10PDI = do.call("rbind", lapply(ResultsPDIM10, "[[", "value")) 
MM10PDISE = do.call("rbind", lapply(ResultsPDIM10, "[[", "se")) 

#COMBINE THE PDI WITH RUBIN'S RULE
PDIcombined <- matrix(ncol = 3, nrow = 1)
colnames(PDIcombined) <- c('PDI', 'LL', 'UL')
PDIcombined <- data.frame(PDIcombined)

PDIcombined$PDI <- mean(MM10PDI)
WithinVar <- mean((MM10PDISE)^2)
BetweenVar <- var(MM10PDI)
PooledVar <- WithinVar + BetweenVar + BetweenVar/NrImp
PDIcombined$PooledSE <- sqrt(PooledVar)

PDIcombined$LL <- PDIcombined$PDI - 1.96*PDIcombined$PooledSE
PDIcombined$UL <- PDIcombined$PDI + 1.96*PDIcombined$PooledSE

PDIcombined


#############
#pairwise AUC

#MLR

#pair 1 benign vs borderline
testdf$bvsb <- testdf$MM1_pred_be / (testdf$MM1_pred_be + testdf$MM1_pred_bo)
AUC.imp(pred = testdf$bvsb, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 0 | testdf$outcome5 == 1,])$Performance

#Pair 2: benign vs stage I
testdf$bvss1 <- testdf$MM1_pred_be / (testdf$MM1_pred_be +testdf$MM1_pred_bs1)
AUC.imp(pred = testdf$bvss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 0 | testdf$outcome5 == 2,])$Performance

#Pair 3: benign vs stage II-IV
testdf$bvss2 <- testdf$MM1_pred_be / (testdf$MM1_pred_be +testdf$MM1_pred_bs2)
AUC.imp(pred = testdf$bvss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 0 | testdf$outcome5 == 3,])$Performance

#Pair 4: benign vs secondary metastatic
testdf$bvsmt <- testdf$MM1_pred_be / (testdf$MM1_pred_be +testdf$MM1_pred_mt)
AUC.imp(pred = testdf$bvss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 0 | testdf$outcome5 == 4,])$Performance

#Pair 5: borderline vs stage I
testdf$bovss1 <- testdf$MM1_pred_bo / (testdf$MM1_pred_bo +testdf$MM1_pred_bs1)
AUC.imp(pred = testdf$bovss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 2,])$Performance

#Pair 6: borderline vs stage II-IV
testdf$bovss2 <- testdf$MM1_pred_bo / (testdf$MM1_pred_bo +testdf$MM1_pred_bs2)
AUC.imp(pred = testdf$bovss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 3,])$Performance

#pair 7: borderline vs secondary metastatic
testdf$bovsmt <- testdf$MM1_pred_bo / (testdf$MM1_pred_bo +testdf$MM1_pred_mt)
AUC.imp(pred = testdf$bovsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 4,])$Performance

#pair 8: stage I vs stage II-IV
testdf$s1vss2 <- testdf$MM1_pred_bs1 / (testdf$MM1_pred_bs1 +testdf$MM1_pred_bs2)
AUC.imp(pred = testdf$s1vss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 3,])$Performance

#pair 9: stage I vs secondary metastatic
testdf$s1vsmt <- testdf$MM1_pred_bs1 / (testdf$MM1_pred_bs1 +testdf$MM1_pred_mt)
AUC.imp(pred = testdf$s1vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 4,])$Performance

#pair 10: stage II-IV vs secondary metastatic
testdf$s2vsmt <- testdf$MM1_pred_bs2 / (testdf$MM1_pred_bs2 +testdf$MM1_pred_mt)
AUC.imp(pred = testdf$s2vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 4,])$Performance


#LINEAR MLR

#pair 1 benign vs borderline
testdf$bvsbm2 <- testdf$MM2_pred_be / (testdf$MM2_pred_be + testdf$MM2_pred_bo)
AUC.imp(pred = testdf$bvsbm2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 0 | testdf$outcome5 == 1,])$Performance

#Pair 2: benign vs stage I
testdf$bvss1m2 <- testdf$MM2_pred_be / (testdf$MM2_pred_be +testdf$MM2_pred_bs1)
AUC.imp(pred = testdf$bvss1m2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 0 | testdf$outcome5 == 2,])$Performance

#Pair 3: benign vs stage II-IV
testdf$bvss2m2 <- testdf$MM2_pred_be / (testdf$MM2_pred_be +testdf$MM2_pred_bs2)
AUC.imp(pred = testdf$bvss2m2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 0 | testdf$outcome5 == 3,])$Performance

#Pair 4: benign vs secondary metastatic
testdf$bvsmtm2 <- testdf$MM2_pred_be / (testdf$MM2_pred_be +testdf$MM2_pred_mt)
AUC.imp(pred = testdf$bvss2m2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 0 | testdf$outcome5 == 4,])$Performance

#Pair 5: borderline vs stage I
testdf$bovss1m2 <- testdf$MM2_pred_bo / (testdf$MM2_pred_bo +testdf$MM2_pred_bs1)
AUC.imp(pred = testdf$bovss1m2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 2,])$Performance

#Pair 6: borderline vs stage II-IV
testdf$bovss2m2 <- testdf$MM2_pred_bo / (testdf$MM2_pred_bo +testdf$MM2_pred_bs2)
AUC.imp(pred = testdf$bovss2m2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 3,])$Performance

#pair 7: borderline vs secondary metastatic
testdf$bovsmtm2 <- testdf$MM2_pred_bo / (testdf$MM2_pred_bo +testdf$MM2_pred_mt)
AUC.imp(pred = testdf$bovsmtm2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
testdf$s1vss2m2 <- testdf$MM2_pred_bs1 / (testdf$MM2_pred_bs1 +testdf$MM2_pred_bs2)
AUC.imp(pred = testdf$s1vss2m2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
testdf$s1vsmtm2 <- testdf$MM2_pred_bs1 / (testdf$MM2_pred_bs1 +testdf$MM2_pred_mt)
AUC.imp(pred = testdf$s1vsmtm2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
testdf$s2vsmtm2 <- testdf$MM2_pred_bs2 / (testdf$MM2_pred_bs2 +testdf$MM2_pred_mt)
AUC.imp(pred = testdf$s2vsmtm2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 4 | testdf$outcome5 == 5,])$Performance


#RIDGE
#pair 1 benign vs borderline
testdf$bvsb <- testdf$MM4_pred_be / (testdf$MM4_pred_be + testdf$MM4_pred_bo)
AUC.imp(pred = testdf$bvsb, outcome = testdf$outcome5, imp =  testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
testdf$bvss1 <- testdf$MM4_pred_be / (testdf$MM4_pred_be +testdf$MM4_pred_bs1)
AUC.imp(pred = testdf$bvss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
testdf$bvss2 <- testdf$MM4_pred_be / (testdf$MM4_pred_be +testdf$MM4_pred_bs2)
AUC.imp(pred = testdf$bvss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
testdf$bvsmt <- testdf$MM4_pred_be / (testdf$MM4_pred_be +testdf$MM4_pred_mt)
AUC.imp(pred = testdf$bvss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
testdf$bovss1 <- testdf$MM4_pred_bo / (testdf$MM4_pred_bo +testdf$MM4_pred_bs1)
AUC.imp(pred = testdf$bovss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
testdf$bovss2 <- testdf$MM4_pred_bo / (testdf$MM4_pred_bo +testdf$MM4_pred_bs2)
AUC.imp(pred = testdf$bovss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
testdf$bovsmt <- testdf$MM4_pred_bo / (testdf$MM4_pred_bo +testdf$MM4_pred_mt)
AUC.imp(pred = testdf$bovsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
testdf$s1vss2 <- testdf$MM4_pred_bs1 / (testdf$MM4_pred_bs1 +testdf$MM4_pred_bs2)
AUC.imp(pred = testdf$s1vss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
testdf$s1vsmt <- testdf$MM4_pred_bs1 / (testdf$MM4_pred_bs1 +testdf$MM4_pred_mt)
AUC.imp(pred = testdf$s1vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
testdf$s2vsmt <- testdf$MM4_pred_bs2 / (testdf$MM4_pred_bs2 +testdf$MM4_pred_mt)
AUC.imp(pred = testdf$s2vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 4 | testdf$outcome5 == 5,])$Performance

#FIRTH
#pair 1 benign vs borderline
testdf$bvsbm5 <- testdf$MM5_pred_be / (testdf$MM5_pred_be + testdf$MM5_pred_bo)
AUC.imp(pred = testdf$bvsbm5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
testdf$bvss1m5 <- testdf$MM5_pred_be / (testdf$MM5_pred_be +testdf$MM5_pred_bs1)
AUC.imp(pred = testdf$bvss1m5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
testdf$bvss2m5 <- testdf$MM5_pred_be / (testdf$MM5_pred_be +testdf$MM5_pred_bs2)
AUC.imp(pred = testdf$bvss2m5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
testdf$bvsmtm5 <- testdf$MM5_pred_be / (testdf$MM5_pred_be +testdf$MM5_pred_mt)
AUC.imp(pred = testdf$bvsmtm5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
testdf$bovss1m5 <- testdf$MM5_pred_bo / (testdf$MM5_pred_bo +testdf$MM5_pred_bs1)
AUC.imp(pred = testdf$bovss1m5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
testdf$bovss2m5 <- testdf$MM5_pred_bo / (testdf$MM5_pred_bo +testdf$MM5_pred_bs2)
AUC.imp(pred = testdf$bovss2m5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
testdf$bovsmtm5 <- testdf$MM5_pred_bo / (testdf$MM5_pred_bo +testdf$MM5_pred_mt)
AUC.imp(pred = testdf$bovsmtm5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
testdf$s1vss2m5 <- testdf$MM5_pred_bs1 / (testdf$MM5_pred_bs1 +testdf$MM5_pred_bs2)
AUC.imp(pred = testdf$s1vss2m5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
testdf$s1vsmtm5 <- testdf$MM5_pred_bs1 / (testdf$MM5_pred_bs1 +testdf$MM5_pred_mt)
AUC.imp(pred = testdf$s1vsmtm5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
testdf$s2vsmtm5 <- testdf$MM5_pred_bs2 / (testdf$MM5_pred_bs2 +testdf$MM5_pred_mt)
AUC.imp(pred = testdf$s2vsmtm5, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 4 | testdf$outcome5 == 5,])$Performance



#rf
#pair 1 benign vs borderline
testdf$bvsb <- testdf$MM7_pred_be / (testdf$MM7_pred_be + testdf$MM7_pred_bo)
AUC.imp(pred = testdf$bvsb, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
testdf$bvss1 <- testdf$MM7_pred_be / (testdf$MM7_pred_be +testdf$MM7_pred_bs1)
AUC.imp(pred = testdf$bvss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
testdf$bvss2 <- testdf$MM7_pred_be / (testdf$MM7_pred_be +testdf$MM7_pred_bs2)
AUC.imp(pred = testdf$bvss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
testdf$bvsmt <- testdf$MM7_pred_be / (testdf$MM7_pred_be +testdf$MM7_pred_mt)
AUC.imp(pred = testdf$bvsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
testdf$bovss1 <- testdf$MM7_pred_bo / (testdf$MM7_pred_bo +testdf$MM7_pred_bs1)
AUC.imp(pred = testdf$bovss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
testdf$bovss2 <- testdf$MM7_pred_bo / (testdf$MM7_pred_bo +testdf$MM7_pred_bs2)
AUC.imp(pred = testdf$bovss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
testdf$bovsmt <- testdf$MM7_pred_bo / (testdf$MM7_pred_bo +testdf$MM7_pred_mt)
AUC.imp(pred = testdf$bovsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
testdf$s1vss2 <- testdf$MM7_pred_bs1 / (testdf$MM7_pred_bs1 +testdf$MM7_pred_bs2)
AUC.imp(pred = testdf$s1vss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
testdf$s1vsmt <- testdf$MM7_pred_bs1 / (testdf$MM7_pred_bs1 +testdf$MM7_pred_mt)
AUC.imp(pred = testdf$s1vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
testdf$s2vsmt <- testdf$MM7_pred_bs2 / (testdf$MM7_pred_bs2 +testdf$MM7_pred_mt)
AUC.imp(pred = testdf$s2vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 4 | testdf$outcome5 == 5,])$Performance

#xgboost
#pair 1 benign vs borderline
testdf$bvsb <- testdf$MM8_pred_be / (testdf$MM8_pred_be + testdf$MM8_pred_bo)
AUC.imp(pred = testdf$bvsb, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
testdf$bvss1 <- testdf$MM8_pred_be / (testdf$MM8_pred_be +testdf$MM8_pred_bs1)
AUC.imp(pred = testdf$bvss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
testdf$bvss2 <- testdf$MM8_pred_be / (testdf$MM8_pred_be +testdf$MM8_pred_bs2)
AUC.imp(pred = testdf$bvss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary mt
testdf$bvsmt <- testdf$MM8_pred_be / (testdf$MM8_pred_be +testdf$MM8_pred_mt)
AUC.imp(pred = testdf$bvsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
testdf$bovss1 <- testdf$MM8_pred_bo / (testdf$MM8_pred_bo +testdf$MM8_pred_bs1)
AUC.imp(pred = testdf$bovss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
testdf$bovss2 <- testdf$MM8_pred_bo / (testdf$MM8_pred_bo +testdf$MM8_pred_bs2)
AUC.imp(pred = testdf$bovss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
testdf$bovsmt <- testdf$MM8_pred_bo / (testdf$MM8_pred_bo +testdf$MM8_pred_mt)
AUC.imp(pred = testdf$bovsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
testdf$s1vss2 <- testdf$MM8_pred_bs1 / (testdf$MM8_pred_bs1 +testdf$MM8_pred_bs2)
AUC.imp(pred = testdf$s1vss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
testdf$s1vsmt <- testdf$MM8_pred_bs1 / (testdf$MM8_pred_bs1 +testdf$MM8_pred_mt)
AUC.imp(pred = testdf$s1vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
testdf$s2vsmt <- testdf$MM8_pred_bs2 / (testdf$MM8_pred_bs2 +testdf$MM8_pred_mt)
AUC.imp(pred = testdf$s2vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 4 | testdf$outcome5 == 5,])$Performance

#nn
#pair 1 benign vs borderline
testdf$bvsb <- testdf$MM9_pred_be / (testdf$MM9_pred_be + testdf$MM9_pred_bo)
AUC.imp(pred = testdf$bvsb, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
testdf$bvss1 <- testdf$MM9_pred_be / (testdf$MM9_pred_be +testdf$MM9_pred_bs1)
AUC.imp(pred = testdf$bvss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
testdf$bvss2 <- testdf$MM9_pred_be / (testdf$MM9_pred_be +testdf$MM9_pred_bs2)
AUC.imp(pred = testdf$bvss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
testdf$bvsmt <- testdf$MM9_pred_be / (testdf$MM9_pred_be +testdf$MM9_pred_mt)
AUC.imp(pred = testdf$bvsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
testdf$bovss1 <- testdf$MM9_pred_bo / (testdf$MM9_pred_bo +testdf$MM9_pred_bs1)
AUC.imp(pred = testdf$bovss1, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
testdf$bovss2 <- testdf$MM9_pred_bo / (testdf$MM9_pred_bo +testdf$MM9_pred_bs2)
AUC.imp(pred = testdf$bovss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
testdf$bovsmt <- testdf$MM9_pred_bo / (testdf$MM9_pred_bo +testdf$MM9_pred_mt)
AUC.imp(pred = testdf$bovsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
testdf$s1vss2 <- testdf$MM9_pred_bs1 / (testdf$MM9_pred_bs1 +testdf$MM9_pred_bs2)
AUC.imp(pred = testdf$s1vss2, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
testdf$s1vsmt <- testdf$MM9_pred_bs1 / (testdf$MM9_pred_bs1 +testdf$MM9_pred_mt)
AUC.imp(pred = testdf$s1vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
testdf$s2vsmt <- testdf$MM9_pred_bs2 / (testdf$MM9_pred_bs2 +testdf$MM9_pred_mt)
AUC.imp(pred = testdf$s2vsmt, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 4 | testdf$outcome5 == 5,])$Performance

#svm
#pair 1 benign vs borderline
testdf$bvsbm10 <- testdf$MM10_pred_be / (testdf$MM10_pred_be + testdf$MM10_pred_bo)
AUC.imp(pred = testdf$bvsbm10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 2,])$Performance

#Pair 2: benign vs stage I
testdf$bvss1m10 <- testdf$MM10_pred_be / (testdf$MM10_pred_be +testdf$MM10_pred_bs1)
AUC.imp(pred = testdf$bvss1m10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 3,])$Performance

#Pair 3: benign vs stage II-IV
testdf$bvss2m10 <- testdf$MM10_pred_be / (testdf$MM10_pred_be +testdf$MM10_pred_bs2)
AUC.imp(pred = testdf$bvss2m10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 4,])$Performance

#Pair 4: benign vs secondary metastatic
testdf$bvsmtm10 <- testdf$MM10_pred_be / (testdf$MM10_pred_be +testdf$MM10_pred_mt)
AUC.imp(pred = testdf$bvss2m10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 1 | testdf$outcome5 == 5,])$Performance

#Pair 5: borderline vs stage I
testdf$bovss1m10 <- testdf$MM10_pred_bo / (testdf$MM10_pred_bo +testdf$MM10_pred_bs1)
AUC.imp(pred = testdf$bovss1m10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 3,])$Performance

#Pair 6: borderline vs stage II-IV
testdf$bovss2m10 <- testdf$MM10_pred_bo / (testdf$MM10_pred_bo +testdf$MM10_pred_bs2)
AUC.imp(pred = testdf$bovss2m10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 4,])$Performance

#pair 7: borderline vs secondary metastatic
testdf$bovsmtm10 <- testdf$MM10_pred_bo / (testdf$MM10_pred_bo +testdf$MM10_pred_mt)
AUC.imp(pred = testdf$bovsmtm10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 2 | testdf$outcome5 == 5,])$Performance

#pair 8: stage I vs stage II-IV
testdf$s1vss2m10 <- testdf$MM10_pred_bs1 / (testdf$MM10_pred_bs1 +testdf$MM10_pred_bs2)
AUC.imp(pred = testdf$s1vss2m10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 4,])$Performance

#pair 9: stage I vs secondary metastatic
testdf$s1vsmtm10 <- testdf$MM10_pred_bs1 / (testdf$MM10_pred_bs1 +testdf$MM10_pred_mt)
AUC.imp(pred = testdf$s1vsmtm10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 3 | testdf$outcome5 == 5,])$Performance

#pair 10: stage II-IV vs secondary metastatic
testdf$s2vsmtm10 <- testdf$MM10_pred_bs2 / (testdf$MM10_pred_bs2 +testdf$MM10_pred_mt)
AUC.imp(pred = testdf$s2vsmtm10, outcome = testdf$outcome5, imp = testdf$imp, data = testdf[testdf$outcome5 == 4 | testdf$outcome5 == 5,])$Performance


#############
# CALIBRATION

####Calibration benign vs malignant

testdf$Center <- factor(testdf$Center, labels=c("Athens", "Milan 1", "Genk", "Leuven", "Katowice", "Malmö", 
                                                "Milan 2", "Monza", "Other", "Pamplona", "Rome", "Cagliari",
                                                "Stockholm", "Trieste"))

IncludedCenters <- unique(testdf$Center)
#selecting line types, colors and width
lty.centers <- c(1,2,3,4,5,6,1,2,3,4,5,6,1,2)
lwd.centers <- c(1,2,1,2,1,2,1,2,1,2,1,2,1,2)

ColorCen <- c("red", "limegreen", "orange", "blue", "maroon", "red", "limegreen", "orange", "blue", "maroon", "red", "limegreen", "orange", "blue")
SampleSize <- c(567, 166, 406, 139, 501, 794, 367, 98, 267, 334, 111, 681, 363, 111)
sqrt(SampleSize)

#LR1
#centre specific calibration
tiff("binary calibration with ca125 mlr.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p=testdf$MM1_mal, y=testdf$outcome1, center=testdf$Center, imputation.id=testdf$imp, patientid=patientid, data=testdf, CL = "CI", CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title="") 
dev.off()
#overall calibration
CM1 <- RE.ValProbImp(p=testdf$MM1_mal, y=testdf$outcome1, center=testdf$Center, imputation.id=testdf$imp, patientid=testdf$patientid, data=testdf)

#linear mlr
#Centre specific
tiff("binary calibration with ca125 linear mlr.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p=testdf$MM2_mal, y=testdf$outcome1, center=testdf$Center, imputation.id=testdf$imp, patientid=patientid, data=testdf, CL = "CI", CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title="") 
dev.off()
#overall
CM2 <- RE.ValProbImp(p=MM2_mal, y=outcome1, center=Center, imputation.id=imp, patientid=patientid, data=testdf)

#RIDGE LR
#centre specific
tiff("binary calibration with ca125 ridge.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p=testdf$MM4_mal, y=testdf$outcome1, center=testdf$Center, imputation.id=testdf$imp, patientid=patientid, data=testdf, CL = "CI", CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title="") 
dev.off()
#overall
CM4 <- RE.ValProbImp(p=MM4_mal, y=outcome1, center=Center, imputation.id=imp, patientid=patientid, data=testdf)

#FIRTH LR
#centre specific
tiff("binary calibration with ca125 firth.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p=testdf$MM5_mal, y=testdf$outcome1, center=testdf$Center, imputation.id=testdf$imp, patientid=patientid, data=testdf, CL = "CI", CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title="") 
dev.off()
#overall
CM5 <- RE.ValProbImp(p=MM5_mal, y=outcome1, center=Center, imputation.id=imp, patientid=patientid, data=testdf)


#RF7
#centre specific
tiff("binary calibration with ca125 rf.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p=testdf$MM7_mal, y=testdf$outcome1, center=testdf$Center, imputation.id=testdf$imp, patientid=patientid, data=testdf, CL = "CI", CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title="") 
dev.off()
#overall
CM7 <- RE.ValProbImp(p=MM7_mal, y=outcome1, center=Center, imputation.id=imp, patientid=patientid, data=testdf)

#XGBOOST
#centre specific
tiff("binary calibration with ca125 xgb.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p=testdf$MM8_mal, y=testdf$outcome1, center=testdf$Center, imputation.id=testdf$imp, patientid=patientid, data=testdf, CL = "CI", CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title="") 
dev.off()
#overall
CM8 <- RE.ValProbImp(p=MM8_mal, y=outcome1, center=Center, imputation.id=imp, patientid=patientid, data=testdf)

#NN
#centre specific
tiff("binary calibration with ca125 nn.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p=testdf$MM9_mal, y=testdf$outcome1, center=testdf$Center, imputation.id=testdf$imp, patientid=patientid, data=testdf, CL = "CI", CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title="") 
dev.off()
#overall
CM9 <- RE.ValProbImp(p=MM9_mal, y=outcome1, center=Center, imputation.id=imp, patientid=patientid, data=testdf)

#SVM
#centre specific
tiff("binary calibration with ca125 svm.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p=testdf$MM10_mal, y=testdf$outcome1, center=testdf$Center, imputation.id=testdf$imp, patientid=patientid, data=testdf, CL = "CI", CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title="")
dev.off()
#overall
CM10 <- RE.ValProbImp(p=MM10_mal, y=outcome1, center=Center, imputation.id=imp, patientid=patientid, data=testdf)


## Predict the outcome for the calibration curve
# LR
OverallCal.M1 = CM1$Plot$y
p.M1 = CM1$Plot$x

# LR2
OverallCal.M2 = CM2$Plot$y
p.M2 = CM2$Plot$x

# ridge
OverallCal.M4 = CM4$Plot$y
p.M4 = CM4$Plot$x

# FIRTH
OverallCal.M5 = CM5$Plot$y
p.M5 = CM5$Plot$x

#RF
OverallCal.M7 = CM7$Plot$y
p.M7 = CM7$Plot$x

# XGBOOST
OverallCal.M8 = CM8$Plot$y
p.M8 = CM8$Plot$x

#NN
OverallCal.M9 = CM9$Plot$y
p.M9 = CM9$Plot$x

# SVM
OverallCal.M10 = CM10$Plot$y
p.M10 = CM10$Plot$x


table <- matrix(ncol = 3, nrow = 8)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'Random forest', 'XGBoost', 'Neural network', 'SVM')
table[1, 2:3] <- c(paste0(format(round(CM1$SumInt$Est, 2), nsmall = 2), " (", format(round(CM1$SumInt$LL, 2), nsmall = 2), "; ", format(round(CM1$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(CM1$SumSlo$Est, 2), nsmall = 2), " (", format(round(CM1$SumSlo$LL, 2), nsmall = 2), "; ", format(round(CM1$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(CM4$SumInt$Est, 2), nsmall = 2), " (", format(round(CM4$SumInt$LL, 2), nsmall = 2), "; ", format(round(CM4$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(CM4$SumSlo$Est, 2), nsmall = 2), " (", format(round(CM4$SumSlo$LL, 2), nsmall = 2), "; ", format(round(CM4$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(CM5$SumInt$Est, 2), nsmall = 2), " (", format(round(CM5$SumInt$LL, 2), nsmall = 2), "; ", format(round(CM5$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(CM5$SumSlo$Est, 2), nsmall = 2), " (", format(round(CM5$SumSlo$LL, 2), nsmall = 2), "; ", format(round(CM5$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(CM2$SumInt$Est, 2), nsmall = 2), " (", format(round(CM2$SumInt$LL, 2), nsmall = 2), "; ", format(round(CM2$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(CM2$SumSlo$Est, 2), nsmall = 2), " (", format(round(CM2$SumSlo$LL, 2), nsmall = 2), "; ", format(round(CM2$SumSlo$UL, 2), nsmall = 2), ")"))
table[5, 2:3] <- c(paste0(format(round(CM7$SumInt$Est, 2), nsmall = 2), " (", format(round(CM7$SumInt$LL, 2), nsmall = 2), "; ", format(round(CM7$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(CM7$SumSlo$Est, 2), nsmall = 2), " (", format(round(CM7$SumSlo$LL, 2), nsmall = 2), "; ", format(round(CM7$SumSlo$UL, 2), nsmall = 2), ")"))
table[6, 2:3] <- c(paste0(format(round(CM8$SumInt$Est, 2), nsmall = 2), " (", format(round(CM8$SumInt$LL, 2), nsmall = 2), "; ", format(round(CM8$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(CM8$SumSlo$Est, 2), nsmall = 2), " (", format(round(CM8$SumSlo$LL, 2), nsmall = 2), "; ", format(round(CM8$SumSlo$UL, 2), nsmall = 2), ")"))
table[7, 2:3] <- c(paste0(format(round(CM9$SumInt$Est, 2), nsmall = 2), " (", format(round(CM9$SumInt$LL, 2), nsmall = 2), "; ", format(round(CM9$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(CM9$SumSlo$Est, 2), nsmall = 2), " (", format(round(CM9$SumSlo$LL, 2), nsmall = 2), "; ", format(round(CM9$SumSlo$UL, 2), nsmall = 2), ")"))
table[8, 2:3] <- c(paste0(format(round(CM10$SumInt$Est, 2), nsmall = 2), " (", format(round(CM10$SumInt$LL, 2), nsmall = 2), "; ", format(round(CM10$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(CM10$SumSlo$Est, 2), nsmall = 2), " (", format(round(CM10$SumSlo$LL, 2), nsmall = 2), "; ", format(round(CM10$SumSlo$UL, 2), nsmall = 2), ")"))

## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Binary calibration per model with ca125.tiff", width = 16, height = 16, units = "cm", res = 300)
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) # cex.lab = 1.5, cex.axis = 1.5
lines(p.M1, OverallCal.M1, lwd = 2, col = "#3399CC")
lines(p.M4, OverallCal.M4, lwd = 2,lty=3, col = "#CCFF66")
lines(p.M5, OverallCal.M5, lwd = 2, lty=5, col = "#FF9933")
lines(p.M2, OverallCal.M2, lwd=2, lty=3, col="#33FFFF")
lines(p.M7, OverallCal.M7, lwd=2,lty=2,  col="#332288")
lines(p.M8, OverallCal.M8, lwd=2, lty=3, col="#00FF66")
lines(p.M9, OverallCal.M9, lwd=2, lty=6, col="#FF33CC")
lines(p.M10, OverallCal.M10, lwd=2, col="#882255")

legend(x = -0.035, y = 1, legend = c( "Ideal", 'MLR', 'Ridge MLR', 'Firth MLR', 'Linear MLR', 'Random forest', 'XGBoost', 'Neural network', 'SVM'),
       col = c("gray50", "#3399CC", "#CCFF66",  "#FF9933","#33FFFF",  "#332288", "#00FF66", "#FF33CC", "#882255"),
       lty = c(2,1,3,5,3,2,3,6,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1)
addtable2plot(x = 0.4, y = -0.02, table = table, display.colnames= TRUE, cex = 0.7) 
dev.off()


################################
#multinomial calibration update#
################################

k <- 5
r <- 1
dfr=2

#MLR

testdfNI = testdf[testdf$imp == 1, ]

pred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

fitpred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

p = as.matrix(testdf[testdf$imp == 1, c('MM1_pred_be', 'MM1_pred_bo', 'MM1_pred_bs1', 'MM1_pred_bs2', 'MM1_pred_mt')])
testdf$Llpbom1 <- logit(testdf$MM1_pred_bo)
testdf$Llps1m1 <- logit(testdf$MM1_pred_bs1)
testdf$Llps2m1 <- logit(testdf$MM1_pred_bs2)
testdf$Llpmtm1 <- logit(testdf$MM1_pred_mt)


LP = as.matrix(testdf[testdf$imp == 1, c('Llpbom1', 'Llps1m1', 'Llps2m1', 'Llpmtm1')])

mean.p <- p
mean.p <- data.frame(mean.p)
mean.p[, 1] <- summaryBy(MM1_pred_be ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 2] <- summaryBy(MM1_pred_bo ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 3] <- summaryBy(MM1_pred_bs1 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 4] <- summaryBy(MM1_pred_bs2 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 5] <- summaryBy(MM1_pred_mt ~ patientid, data = testdf, FUN = mean)[,2]


for(i in 1:10){
  testImpLR1 <- testdf[testdf$imp == i, ]
  
  
  p = as.matrix(testImpLR1[, c('MM1_pred_be', 'MM1_pred_bo', 'MM1_pred_bs1', 'MM1_pred_bs2', 'MM1_pred_mt')])
  LP = as.matrix(testImpLR1[, c('Llpbom1', 'Llps1m1', 'Llps2m1', 'Llpmtm1')])
  
  lps <- split(LP,col(LP))
  for(j in 1:(k-1)){assign(paste("lp", j, sep = ""),unlist(lps[[j]]))}
  
  outcome = testImpLR1$outcome5
  r=1
  fitp<-vglm(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)
  
  pred1[, i] <- p[,1]
  pred2[, i] <- p[,2]
  pred3[, i] <- p[,3]
  pred4[, i] <- p[,4]
  pred5[, i] <- p[,5]
  
  fitpred1[, i] <- fitted(fitp)[,1]
  fitpred2[, i] <- fitted(fitp)[,2]
  fitpred3[, i] <- fitted(fitp)[,3]
  fitpred4[, i] <- fitted(fitp)[,4]
  fitpred5[, i] <- fitted(fitp)[,5]
  
  
}

LP.predict <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(LP.predict) <- c('1', '2', '3', '4', '5')
LP.predict[, 1] <- rowMeans(pred1)
LP.predict[, 2] <- rowMeans(pred2)
LP.predict[, 3] <- rowMeans(pred3)
LP.predict[, 4] <- rowMeans(pred4)
LP.predict[, 5] <- rowMeans(pred5)

FIT <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(FIT) <- c('1', '2', '3', '4', '5')
FIT[, 1] <- rowMeans(fitpred1)
FIT[, 2] <- rowMeans(fitpred2)
FIT[, 3] <- rowMeans(fitpred3)
FIT[, 4] <- rowMeans(fitpred4)
FIT[, 5] <- rowMeans(fitpred5)

probs <- split(LP.predict,col(LP.predict))
FIT <- data.frame(FIT)

dfmlr <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfmlr) <- c('p1', 'fitted', 'Tumour')
dfmlr <- data.frame(dfmlr)
dfmlr$p1 <- unlist(probs)
dfmlr$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
dfmlr$fitted <- stack(FIT[,1:5])
dfmlr$Tumour <- factor(dfmlr$Tumour, levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))

p1 <-ggplot(dfmlr)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability of malignancy")+
  ylab("Observed proportion of malignancy")+geom_text(x=0.1, y=1, label="MLR")
p1

#Ridge MLR

pred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

fitpred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

p = as.matrix(testdf[testdf$imp == 1, c('MM4_pred_be', 'MM4_pred_bo', 'MM4_pred_bs1', 'MM4_pred_bs2', 'MM4_pred_mt')])
testdf$Llpbom4 <- car::logit(testdf$MM4_pred_bo)
testdf$Llps1m4 <- car::logit(testdf$MM4_pred_bs1)
testdf$Llps2m4 <- car::logit(testdf$MM4_pred_bs2)
testdf$Llpmtm4 <- car::logit(testdf$MM4_pred_mt)


LP = as.matrix(testdf[testdf$imp == 1, c('Llpbom4', 'Llps1m4', 'Llps2m4', 'Llpmtm4')])

mean.p <- p
mean.p <- data.frame(mean.p)
mean.p[, 1] <- summaryBy(MM4_pred_be ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 2] <- summaryBy(MM4_pred_bo ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 3] <- summaryBy(MM4_pred_bs1 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 4] <- summaryBy(MM4_pred_bs2 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 5] <- summaryBy(MM4_pred_mt ~ patientid, data = testdf, FUN = mean)[,2]

for(i in 1:10){
  testImpLR4 <- testdf[testdf$imp == i, ]
  
  
  p = as.matrix(testImpLR4[, c('MM4_pred_be', 'MM4_pred_bo', 'MM4_pred_bs1', 'MM4_pred_bs2', 'MM4_pred_mt')])
  LP = as.matrix(testImpLR4[, c('Llpbom4', 'Llps1m4', 'Llps2m4', 'Llpmtm4')])
  
  outcome = testImpLR4$outcome5
  r=1
  
  
  lps <- split(LP,col(LP))
  for(j in 1:(k-1)){assign(paste("lp", j, sep = ""),unlist(lps[[j]]))}
  
  fitp<-vglm(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)
  
  pred1[, i] <- p[,1]
  pred2[, i] <- p[,2]
  pred3[, i] <- p[,3]
  pred4[, i] <- p[,4]
  pred5[, i] <- p[,5]
  
  fitpred1[, i] <- fitted(fitp)[,1]
  fitpred2[, i] <- fitted(fitp)[,2]
  fitpred3[, i] <- fitted(fitp)[,3]
  fitpred4[, i] <- fitted(fitp)[,4]
  fitpred5[, i] <- fitted(fitp)[,5]
  
  
}

LP.predict <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(LP.predict) <- c('1', '2', '3', '4', '5')
LP.predict[, 1] <- rowMeans(pred1)
LP.predict[, 2] <- rowMeans(pred2)
LP.predict[, 3] <- rowMeans(pred3)
LP.predict[, 4] <- rowMeans(pred4)
LP.predict[, 5] <- rowMeans(pred5)

FIT <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(FIT) <- c('1', '2', '3', '4', '5')
FIT[, 1] <- rowMeans(fitpred1)
FIT[, 2] <- rowMeans(fitpred2)
FIT[, 3] <- rowMeans(fitpred3)
FIT[, 4] <- rowMeans(fitpred4)
FIT[, 5] <- rowMeans(fitpred5)

probs <- split(LP.predict,col(LP.predict)) 
FIT <- data.frame(FIT)

dfmlr <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfmlr) <- c('p1', 'fitted', 'Tumour')
dfmlr <- data.frame(dfmlr)
dfmlr$p1 <- unlist(probs)
dfmlr$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
dfmlr$fitted <- stack(FIT[,1:5])
dfmlr$Tumour <- factor(dfmlr$Tumour, levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))

p2 <-ggplot(dfmlr)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta"))+geom_text(x=0.1, y=1, label="Ridge MLR")
p2

#firth

pred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

fitpred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

p = as.matrix(testdf[testdf$imp == 1, c('MM5_pred_be', 'MM5_pred_bo', 'MM5_pred_bs1', 'MM5_pred_bs2', 'MM5_pred_mt')])
testdf$Llpbom5 <- car::logit(testdf$MM5_pred_bo)
testdf$Llps1m5 <- car::logit(testdf$MM5_pred_bs1)
testdf$Llps2m5 <- car::logit(testdf$MM5_pred_bs2)
testdf$Llpmtm5 <- car::logit(testdf$MM5_pred_mt)


LP = as.matrix(testdf[testdf$imp == 1, c('Llpbom5', 'Llps1m5', 'Llps2m5', 'Llpmtm5')])

mean.p <- p
mean.p <- data.frame(mean.p)
mean.p[, 1] <- summaryBy(MM5_pred_be ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 2] <- summaryBy(MM5_pred_bo ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 3] <- summaryBy(MM5_pred_bs1 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 4] <- summaryBy(MM5_pred_bs2 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 5] <- summaryBy(MM5_pred_mt ~ patientid, data = testdf, FUN = mean)[,2]

for(i in 1:10){
  testImpLR5 <- testdf[testdf$imp == i, ]
  
  
  
  
  p = as.matrix(testImpLR5[, c('MM5_pred_be', 'MM5_pred_bo', 'MM5_pred_bs1', 'MM5_pred_bs2', 'MM5_pred_mt')])
  LP = as.matrix(testImpLR5[, c('Llpbom5', 'Llps1m5', 'Llps2m5', 'Llpmtm5')])
  
  outcome = testImpLR5$outcome5
  r=1
  
  lps <- split(LP,col(LP))
  for(j in 1:(k-1)){assign(paste("lp", j, sep = ""),unlist(lps[[j]]))}
  
  fitp<-vglm(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)
  
  #fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
  
  pred1[, i] <- p[,1]
  pred2[, i] <- p[,2]
  pred3[, i] <- p[,3]
  pred4[, i] <- p[,4]
  pred5[, i] <- p[,5]
  
  fitpred1[, i] <- fitted(fitp)[,1]
  fitpred2[, i] <- fitted(fitp)[,2]
  fitpred3[, i] <- fitted(fitp)[,3]
  fitpred4[, i] <- fitted(fitp)[,4]
  fitpred5[, i] <- fitted(fitp)[,5]
  
  
}

LP.predict <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(LP.predict) <- c('1', '2', '3', '4', '5')
LP.predict[, 1] <- rowMeans(pred1)
LP.predict[, 2] <- rowMeans(pred2)
LP.predict[, 3] <- rowMeans(pred3)
LP.predict[, 4] <- rowMeans(pred4)
LP.predict[, 5] <- rowMeans(pred5)

FIT <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(FIT) <- c('1', '2', '3', '4', '5')
FIT[, 1] <- rowMeans(fitpred1)
FIT[, 2] <- rowMeans(fitpred2)
FIT[, 3] <- rowMeans(fitpred3)
FIT[, 4] <- rowMeans(fitpred4)
FIT[, 5] <- rowMeans(fitpred5)

probs <- split(LP.predict,col(LP.predict))
FIT <- data.frame(FIT)

dfmlr <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfmlr) <- c('p1', 'fitted', 'Tumour')
dfmlr <- data.frame(dfmlr)
dfmlr$p1 <- unlist(probs)
dfmlr$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
dfmlr$fitted <- stack(FIT[,1:5])
dfmlr$Tumour <- factor(dfmlr$Tumour, levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))

p3 <-ggplot(dfmlr)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability of malignancy")+
  ylab("Observed proportion of malignancy")+geom_text(x=0.1, y=1, label="Firth MLR")
p3

#linear MLR

pred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

fitpred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

p = as.matrix(testdf[testdf$imp == 1, c('MM2_pred_be', 'MM2_pred_bo', 'MM2_pred_bs1', 'MM2_pred_bs2', 'MM2_pred_mt')])
testdf$Llpbom2 <- logit(testdf$MM2_pred_bo)
testdf$Llps1m2 <- logit(testdf$MM2_pred_bs1)
testdf$Llps2m2 <- logit(testdf$MM2_pred_bs2)
testdf$Llpmtm2 <- logit(testdf$MM2_pred_mt)


LP = as.matrix(testdf[testdf$imp == 1, c('Llpbom2', 'Llps1m2', 'Llps2m2', 'Llpmtm2')])

mean.p <- p
mean.p <- data.frame(mean.p)
mean.p[, 1] <- summaryBy(MM2_pred_be ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 2] <- summaryBy(MM2_pred_bo ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 3] <- summaryBy(MM2_pred_bs1 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 4] <- summaryBy(MM2_pred_bs2 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 5] <- summaryBy(MM2_pred_mt ~ patientid, data = testdf, FUN = mean)[,2]

for(i in 1:10){
  testImpLR2 <- testdf[testdf$imp == i, ]
  
  
  p = as.matrix(testImpLR2[, c('MM2_pred_be', 'MM2_pred_bo', 'MM2_pred_bs1', 'MM2_pred_bs2', 'MM2_pred_mt')])
  LP = as.matrix(testImpLR2[, c('Llpbom2', 'Llps1m2', 'Llps2m2', 'Llpmtm2')])
  
  outcome = testImpLR2$outcome5
  r=1
  #fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
  
  lps <- split(LP,col(LP))
  for(j in 1:(k-1)){assign(paste("lp", j, sep = ""),unlist(lps[[j]]))}
  
  
  fitp<-vglm(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)
  
  pred1[, i] <- p[,1]
  pred2[, i] <- p[,2]
  pred3[, i] <- p[,3]
  pred4[, i] <- p[,4]
  pred5[, i] <- p[,5]
  
  fitpred1[, i] <- fitted(fitp)[,1]
  fitpred2[, i] <- fitted(fitp)[,2]
  fitpred3[, i] <- fitted(fitp)[,3]
  fitpred4[, i] <- fitted(fitp)[,4]
  fitpred5[, i] <- fitted(fitp)[,5]
  
  
}

LP.predict <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(LP.predict) <- c('1', '2', '3', '4', '5')
LP.predict[, 1] <- rowMeans(pred1)
LP.predict[, 2] <- rowMeans(pred2)
LP.predict[, 3] <- rowMeans(pred3)
LP.predict[, 4] <- rowMeans(pred4)
LP.predict[, 5] <- rowMeans(pred5)

FIT <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(FIT) <- c('1', '2', '3', '4', '5')
FIT[, 1] <- rowMeans(fitpred1)
FIT[, 2] <- rowMeans(fitpred2)
FIT[, 3] <- rowMeans(fitpred3)
FIT[, 4] <- rowMeans(fitpred4)
FIT[, 5] <- rowMeans(fitpred5)

probs <- split(LP.predict,col(LP.predict)) 

FIT <- data.frame(FIT)

dfmlr <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfmlr) <- c('p1', 'fitted', 'Tumour')
dfmlr <- data.frame(dfmlr)
dfmlr$p1 <- unlist(probs)
dfmlr$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
dfmlr$fitted <- stack(FIT[,1:5])
dfmlr$Tumour <- factor(dfmlr$Tumour, levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))

p4 <-ggplot(dfmlr)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability ")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="Linear MLR")
p4

.#rf

pred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

fitpred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

p = as.matrix(testdf[testdf$imp == 1, c('MM7_pred_be', 'MM7_pred_bo', 'MM7_pred_bs1', 'MM7_pred_bs2', 'MM7_pred_mt')])
testdf$Llpbom7 <- car::logit(testdf$MM7_pred_bo)
testdf$Llps1m7 <- car::logit(testdf$MM7_pred_bs1)
testdf$Llps2m7 <- car::logit(testdf$MM7_pred_bs2)
testdf$Llpmtm7 <- car::logit(testdf$MM7_pred_mt)


LP = as.matrix(testdf[testdf$imp == 1, c('Llpbom7', 'Llps1m7', 'Llps2m7', 'Llpmtm7')])

mean.p <- p
mean.p <- data.frame(mean.p)
mean.p[, 1] <- summaryBy(MM7_pred_be ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 2] <- summaryBy(MM7_pred_bo ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 3] <- summaryBy(MM7_pred_bs1 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 4] <- summaryBy(MM7_pred_bs2 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 5] <- summaryBy(MM7_pred_mt ~ patientid, data = testdf, FUN = mean)[,2]

for(i in 1:10){
  testImpLR7 <- testdf[testdf$imp == i, ]
  
  
  
  
  p = as.matrix(testImpLR7[, c('MM7_pred_be', 'MM7_pred_bo', 'MM7_pred_bs1', 'MM7_pred_bs2', 'MM7_pred_mt')])
  LP = as.matrix(testImpLR7[, c('Llpbom7', 'Llps1m7', 'Llps2m7', 'Llpmtm7')])
  
  outcome = testImpLR7$outcome5
  r=1
  fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
  
  pred1[, i] <- p[,1]
  pred2[, i] <- p[,2]
  pred3[, i] <- p[,3]
  pred4[, i] <- p[,4]
  pred5[, i] <- p[,5]
  
  fitpred1[, i] <- fitted(fitp)[,1]
  fitpred2[, i] <- fitted(fitp)[,2]
  fitpred3[, i] <- fitted(fitp)[,3]
  fitpred4[, i] <- fitted(fitp)[,4]
  fitpred5[, i] <- fitted(fitp)[,5]
  
  
}

LP.predict <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(LP.predict) <- c('1', '2', '3', '4', '5')
LP.predict[, 1] <- rowMeans(pred1)
LP.predict[, 2] <- rowMeans(pred2)
LP.predict[, 3] <- rowMeans(pred3)
LP.predict[, 4] <- rowMeans(pred4)
LP.predict[, 5] <- rowMeans(pred5)

FIT <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(FIT) <- c('1', '2', '3', '4', '5')
FIT[, 1] <- rowMeans(fitpred1)
FIT[, 2] <- rowMeans(fitpred2)
FIT[, 3] <- rowMeans(fitpred3)
FIT[, 4] <- rowMeans(fitpred4)
FIT[, 5] <- rowMeans(fitpred5)

probs <- split(LP.predict,col(LP.predict)) 

FIT <- data.frame(FIT)

dfmlr <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfmlr) <- c('p1', 'fitted', 'Tumour')
dfmlr <- data.frame(dfmlr)
dfmlr$p1 <- unlist(probs)
dfmlr$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
dfmlr$fitted <- stack(FIT[,1:5])
dfmlr$Tumour <- factor(dfmlr$Tumour, levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))

p5 <-ggplot(dfmlr)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability of malignancy")+
  ylab("Observed proportion of malignancy")+geom_text(x=0.1, y=1, label="RF")
p5

#xgboost

pred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

fitpred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

p = as.matrix(testdf[testdf$imp == 1, c('MM8_pred_be', 'MM8_pred_bo', 'MM8_pred_bs1', 'MM8_pred_bs2', 'MM8_pred_mt')])
testdf$Llpbom8 <- car::logit(testdf$MM8_pred_bo)
testdf$Llps1m8 <- car::logit(testdf$MM8_pred_bs1)
testdf$Llps2m8 <- car::logit(testdf$MM8_pred_bs2)
testdf$Llpmtm8 <- car::logit(testdf$MM8_pred_mt)


LP = as.matrix(testdf[testdf$imp == 1, c('Llpbom8', 'Llps1m8', 'Llps2m8', 'Llpmtm8')])

mean.p <- p
mean.p <- data.frame(mean.p)
mean.p[, 1] <- summaryBy(MM8_pred_be ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 2] <- summaryBy(MM8_pred_bo ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 3] <- summaryBy(MM8_pred_bs1 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 4] <- summaryBy(MM8_pred_bs2 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 5] <- summaryBy(MM8_pred_mt ~ patientid, data = testdf, FUN = mean)[,2]

for(i in 1:10){
  testImpLR8 <- testdf[testdf$imp == i, ]
  
  
  
  
  p = as.matrix(testImpLR8[, c('MM8_pred_be', 'MM8_pred_bo', 'MM8_pred_bs1', 'MM8_pred_bs2', 'MM8_pred_mt')])
  LP = as.matrix(testImpLR8[, c('Llpbom8', 'Llps1m8', 'Llps2m8', 'Llpmtm8')])
  
  outcome = testImpLR8$outcome5
  r=1
  fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
  
  pred1[, i] <- p[,1]
  pred2[, i] <- p[,2]
  pred3[, i] <- p[,3]
  pred4[, i] <- p[,4]
  pred5[, i] <- p[,5]
  
  fitpred1[, i] <- fitted(fitp)[,1]
  fitpred2[, i] <- fitted(fitp)[,2]
  fitpred3[, i] <- fitted(fitp)[,3]
  fitpred4[, i] <- fitted(fitp)[,4]
  fitpred5[, i] <- fitted(fitp)[,5]
  
  
}

LP.predict <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(LP.predict) <- c('1', '2', '3', '4', '5')
LP.predict[, 1] <- rowMeans(pred1)
LP.predict[, 2] <- rowMeans(pred2)
LP.predict[, 3] <- rowMeans(pred3)
LP.predict[, 4] <- rowMeans(pred4)
LP.predict[, 5] <- rowMeans(pred5)

FIT <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(FIT) <- c('1', '2', '3', '4', '5')
FIT[, 1] <- rowMeans(fitpred1)
FIT[, 2] <- rowMeans(fitpred2)
FIT[, 3] <- rowMeans(fitpred3)
FIT[, 4] <- rowMeans(fitpred4)
FIT[, 5] <- rowMeans(fitpred5)

probs <- split(LP.predict,col(LP.predict)) 
FIT <- data.frame(FIT)

dfmlr <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfmlr) <- c('p1', 'fitted', 'Tumour')
dfmlr <- data.frame(dfmlr)
dfmlr$p1 <- unlist(probs)
dfmlr$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
dfmlr$fitted <- stack(FIT[,1:5])
dfmlr$Tumour <- factor(dfmlr$Tumour, levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))

p6 <-ggplot(dfmlr)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability ")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="XGBoost")
p6
#nn

pred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

fitpred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

p = as.matrix(testdf[testdf$imp == 1, c('MM9_pred_be', 'MM9_pred_bo', 'MM9_pred_bs1', 'MM9_pred_bs2', 'MM9_pred_mt')])
testdf$Llpbom9 <- car::logit(testdf$MM9_pred_bo)
testdf$Llps1m9 <- car::logit(testdf$MM9_pred_bs1)
testdf$Llps2m9 <- car::logit(testdf$MM9_pred_bs2)
testdf$Llpmtm9 <- car::logit(testdf$MM9_pred_mt)


LP = as.matrix(testdf[testdf$imp == 1, c('Llpbom9', 'Llps1m9', 'Llps2m9', 'Llpmtm9')])

mean.p <- p
mean.p <- data.frame(mean.p)
mean.p[, 1] <- summaryBy(MM9_pred_be ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 2] <- summaryBy(MM9_pred_bo ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 3] <- summaryBy(MM9_pred_bs1 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 4] <- summaryBy(MM9_pred_bs2 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 5] <- summaryBy(MM9_pred_mt ~ patientid, data = testdf, FUN = mean)[,2]

for(i in 1:10){
  testImpLR9 <- testdf[testdf$imp == i, ]
  
  
  
  
  p = as.matrix(testImpLR9[, c('MM9_pred_be', 'MM9_pred_bo', 'MM9_pred_bs1', 'MM9_pred_bs2', 'MM9_pred_mt')])
  LP = as.matrix(testImpLR9[, c('Llpbom9', 'Llps1m9', 'Llps2m9', 'Llpmtm9')])
  
  outcome = testImpLR9$outcome5
  r=1
  fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
  
  pred1[, i] <- p[,1]
  pred2[, i] <- p[,2]
  pred3[, i] <- p[,3]
  pred4[, i] <- p[,4]
  pred5[, i] <- p[,5]
  
  fitpred1[, i] <- fitted(fitp)[,1]
  fitpred2[, i] <- fitted(fitp)[,2]
  fitpred3[, i] <- fitted(fitp)[,3]
  fitpred4[, i] <- fitted(fitp)[,4]
  fitpred5[, i] <- fitted(fitp)[,5]
  
  
}

LP.predict <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(LP.predict) <- c('1', '2', '3', '4', '5')
LP.predict[, 1] <- rowMeans(pred1)
LP.predict[, 2] <- rowMeans(pred2)
LP.predict[, 3] <- rowMeans(pred3)
LP.predict[, 4] <- rowMeans(pred4)
LP.predict[, 5] <- rowMeans(pred5)

FIT <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(FIT) <- c('1', '2', '3', '4', '5')
FIT[, 1] <- rowMeans(fitpred1)
FIT[, 2] <- rowMeans(fitpred2)
FIT[, 3] <- rowMeans(fitpred3)
FIT[, 4] <- rowMeans(fitpred4)
FIT[, 5] <- rowMeans(fitpred5)

probs <- split(LP.predict,col(LP.predict)) 
FIT <- data.frame(FIT)

dfmlr <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfmlr) <- c('p1', 'fitted', 'Tumour')
dfmlr <- data.frame(dfmlr)
dfmlr$p1 <- unlist(probs)
dfmlr$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
dfmlr$fitted <- stack(FIT[,1:5])
dfmlr$Tumour <- factor(dfmlr$Tumour, levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))

p7 <-ggplot(dfmlr)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_rect(color = "black",fill = NA))+scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta")) +
  xlab("Estimated probability of malignancy")+
  ylab("Observed proportion of malignancy")+geom_text(x=0.1, y=1, label="NN")
p7

#svm
coefp <- matrix(nrow = 10, ncol = 20)
colnames(coefp) <- c('(Intercept):1', '(Intercept):2', '(Intercept):3', '(Intercept):4',
                     'LPlp2w:1', 'LPlp2w:2', 'LPlp2w:3', 'LPlp2w:4',
                     'LPlp3w:1', 'LPlp3w:2', 'LPlp3w:3', 'LPlp3w:4',
                     'LPlp4w:1', 'LPlp4w:2', 'LPlp4w:3', 'LPlp4w:4',
                     'LPlp5w:1', 'LPlp5w:2', 'LPlp5w:3', 'LPlp5w:4')
pred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
pred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

fitpred1 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred2 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred3 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred4 <- matrix(nrow = nrow(testdfNI), ncol = 10)
fitpred5 <- matrix(nrow = nrow(testdfNI), ncol = 10)

p = as.matrix(testdf[testdf$imp == 1, c('MM10_pred_be', 'MM10_pred_bo', 'MM10_pred_bs1', 'MM10_pred_bs2', 'MM10_pred_mt')])
testdf$Llpbom10 <- car::logit(testdf$MM10_pred_bo)
testdf$Llps1m10 <- car::logit(testdf$MM10_pred_bs1)
testdf$Llps2m10 <- car::logit(testdf$MM10_pred_bs2)
testdf$Llpmtm10 <- car::logit(testdf$MM10_pred_mt)


LP = as.matrix(testdf[testdf$imp == 1, c('Llpbom10', 'Llps1m10', 'Llps2m10', 'Llpmtm10')])

mean.p <- p
mean.p <- data.frame(mean.p)
mean.p[, 1] <- summaryBy(MM10_pred_be ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 2] <- summaryBy(MM10_pred_bo ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 3] <- summaryBy(MM10_pred_bs1 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 4] <- summaryBy(MM10_pred_bs2 ~ patientid, data = testdf, FUN = mean)[,2]
mean.p[, 5] <- summaryBy(MM10_pred_mt ~ patientid, data = testdf, FUN = mean)[,2]

for(i in 1:10){
  testImpLR10 <- testdf[testdf$imp == i, ]
  
  
  p = as.matrix(testImpLR10[, c('MM10_pred_be', 'MM10_pred_bo', 'MM10_pred_bs1', 'MM10_pred_bs2', 'MM10_pred_mt')])
  LP = as.matrix(testImpLR10[, c('Llpbom10', 'Llps1m10', 'Llps2m10', 'Llpmtm10')])
  
  outcome = testImpLR10$outcome5
  r=1
  fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
  
  pred1[, i] <- p[,1]
  pred2[, i] <- p[,2]
  pred3[, i] <- p[,3]
  pred4[, i] <- p[,4]
  pred5[, i] <- p[,5]
  
  fitpred1[, i] <- fitted(fitp)[,1]
  fitpred2[, i] <- fitted(fitp)[,2]
  fitpred3[, i] <- fitted(fitp)[,3]
  fitpred4[, i] <- fitted(fitp)[,4]
  fitpred5[, i] <- fitted(fitp)[,5]
  
  
}

LP.predict <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(LP.predict) <- c('1', '2', '3', '4', '5')
LP.predict[, 1] <- rowMeans(pred1)
LP.predict[, 2] <- rowMeans(pred2)
LP.predict[, 3] <- rowMeans(pred3)
LP.predict[, 4] <- rowMeans(pred4)
LP.predict[, 5] <- rowMeans(pred5)

FIT <- matrix(nrow = nrow(testdfNI), ncol = 5)
colnames(FIT) <- c('1', '2', '3', '4', '5')
FIT[, 1] <- rowMeans(fitpred1)
FIT[, 2] <- rowMeans(fitpred2)
FIT[, 3] <- rowMeans(fitpred3)
FIT[, 4] <- rowMeans(fitpred4)
FIT[, 5] <- rowMeans(fitpred5)

probs <- split(LP.predict,col(LP.predict)) 
FIT <- data.frame(FIT)

dfmlr <- matrix(ncol = 3, nrow = 2489*5)
colnames(dfmlr) <- c('p1', 'fitted', 'Tumour')
dfmlr <- data.frame(dfmlr)
dfmlr$p1 <- unlist(probs)
dfmlr$Tumour <- rep(c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'), each=2489)
dfmlr$fitted <- stack(FIT[,1:5])
dfmlr$Tumour <- factor(dfmlr$Tumour, levels=c('Benign', 'Borderline', 'Stage I invasive', 'Stage II-IV invasive', 'Secondary metastatic'))

p8 <-ggplot(dfmlr)+geom_spline(aes(x=p1, y=fitted$values, colour=Tumour, group=Tumour),df=2, size=1.3)+geom_abline(intercept =0,slope = 1,size=1.3, colour="gray",linetype=2)+
  xlab("Estimated probability of malignancy")+
  ylab("Observed proportion")+geom_text(x=0.1, y=1, label="SVM")+
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),panel.grid.minor = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),panel.border = element_rect(color = "black",fill = NA))+
  scale_color_manual(values=c("red", "limegreen", "darkblue", "cyan", "magenta"))
p8

#make plot
tiff("multinomial calibration with CA125.tiff", width = 22, height = 31, units = "cm", res = 300)
figure <- ggarrange(p1, p2, p3, p4, p5,p6,p7,p8, ncol=2, nrow=4, common.legend=TRUE, legend="bottom", align="hv")
figure
dev.off()

##################################
## CLINICAL UTILITY NET BENEFIT ##
##################################

CULR1 <- DataWinBugs.imp(pred=MM1_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)
CULR4 <- DataWinBugs.imp(pred=MM4_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)
CULR5 <- DataWinBugs.imp(pred=MM5_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)
CULR2 <- DataWinBugs.imp(pred=MM2_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)
CULR7 <- DataWinBugs.imp(pred=MM7_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)
CULR8 <- DataWinBugs.imp(pred=MM8_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)
CULR9 <- DataWinBugs.imp(pred=MM9_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)
CULR10 <- DataWinBugs.imp(pred=MM10_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)


#Plot
tiff("NB with CA125.tiff", width = 16, height = 16, units = "cm", res = 300)

plot(nb$threshold,nb$`NB treatall`, type = "l", col = "gray50", lty=1, xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.01, 0.35), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(nb$threshold, nb$`NB LR1`, type='l', col = "#3399CC", lty = 2, lwd = 2)
points(nb$threshold,nb$`NB RIDGE`, type = 'l', col = "#CCFF66", lty = 3, lwd = 3)
points(nb$threshold, nb$`NB FIRTH`, type = 'l', col = "#FF9933", lty=4,lwd = 1)
points(nb$threshold, nb$`NB LR2`, type = 'l', col = "#33FFFF", lty=5,lwd = 2)
points(nb$threshold, nb$`NB RF`, type = 'l', col = "#332288", lty = 6, lwd = 3)
points(nb$threshold, nb$`NB XGB`, type = 'l', col = "#00FF66", lty=1,lwd = 1)
points(nb$threshold, nb$`NB NN`, type = 'l', col = "#FF33CC", lty=2,lwd = 2)
points(nb$threshold, nb$`NB SVM`, type='l', col = "#882255", lty = 3, lwd = 3)
abline(a = 0, b = 0, lty = 4, lwd = 1, col = "gray50")
legend(x=0.05, y=0.15,legend = c("Treat all","MLR", "Ridge MLR", "Firth MLR", "Linear MLR", "Random forest", "XGBoost", "Neural network", "Support vector machine","Treat none"), ncol = 1, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("gray50", "#3399CC", "#CCFF66",  "#FF9933", "#33FFFF", '#332288', '#00FF66','#FF33CC',"#882255", 'gray50'), 
       lty = c(1,2,3,4,5,6,1,2,3,4), lwd = c(1,2,3,1,2,3,1,2,3,1), cex = 0.8, bty = "n")
dev.off()


## NET BENEFIT FOR MULTINOMIAL OUTCOMES

#Benign vs the rest
testdf$benign <- as.factor(ifelse(testdf$outcome5 == 1, 1, 0))
#MLR MM1pred

m1benign <- NB.imp(pred=MM1_pred_be, outcome=benign, data=testdf, imp=imp)$Results
m1benign

#ridge MLR
m4benign <- NB.imp(pred=MM4_pred_be, outcome=benign, data=testdf, imp=imp)$Results
m4benign
#firth MLR
m5benign <- NB.imp(pred=MM5_pred_be, outcome=benign, data=testdf, imp=imp)$Results
m5benign
#linear MLR
m2benign <- NB.imp(pred=MM2_pred_be, outcome=benign, data=testdf, imp=imp)$Results
m2benign
#RF
m7benign <- NB.imp(pred=MM7_pred_be, outcome=benign, data=testdf, imp=imp)$Results
m7benign
#XGBoost
m8benign <- NB.imp(pred=MM8_pred_be, outcome=benign, data=testdf, imp=imp)$Results
m8benign
#NN
m9benign <- NB.imp(pred=MM9_pred_be, outcome=benign, data=testdf, imp=imp)$Results
m9benign
#SVM
m10benign <- NB.imp(pred=MM10_pred_be, outcome=benign, data=testdf, imp=imp)$Results
m10benign

benignnb <- data.frame(m1benign$CutOff,m1benign$NB.mean, m2benign$NB.mean, m4benign$NB.mean, m5benign$NB.mean, m7benign$NB.mean, m8benign$NB.mean, m9benign$NB.mean, m10benign$NB.mean, m1benign$NBTA.mean)

#plot
tiff("netbenefit benign.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(benignnb$m1benign.CutOff,benignnb$m1benign.NBTA.mean, type = "l", col = "#440154FF", lty=1, xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-1, 0), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(benignnb$m1benign.CutOff, benignnb$m1benign.NB.mean, type='l', col = "#482677FF", lty = 2, lwd = 2)
points(benignnb$m1benign.CutOff, benignnb$m4benign.NB.mean, type = 'l', col = "#404788FF", lty = 3, lwd = 3)
points(benignnb$m1benign.CutOff, benignnb$m5benign.NB.mean, type = 'l', col = "#33638DFF", lty = 4, lwd = 1)
points(benignnb$m1benign.CutOff, benignnb$m2benign.NB.mean, type = 'l', col = "#238A8DFF", lty=5,lwd = 2)
points(benignnb$m1benign.CutOff, benignnb$m7benign.NB.mean, type = 'l', col = "#20A387FF", lty=6,lwd = 3)
points(benignnb$m1benign.CutOff, benignnb$m8benign.NB.mean, type = 'l', col = "#55C667FF", lty = 1, lwd = 1)
points(benignnb$m1benign.CutOff, benignnb$m9benign.NB.mean, type = 'l', col = "#95D840FF", lty = 2, lwd = 2)
points(benignnb$m1benign.CutOff, benignnb$m10benign.NB.mean, type = 'l', col = "#FDE725FF", lty=4,lwd = 1)
abline(a = 0, b = 0, lty = 5, lwd = 1, col = "gray50")
legend("bottomleft", legend = c( "Treat all", "MLR", "Ridge MLR", "Firth MLR", "Linear MLR", "Random forest", "XGBoost", "Neural network", "Support vector machine", "Treat none"), ncol = 1, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("#440154FF", "#482677FF", "#404788FF",  "#33638DFF", "#238A8DFF", '#20A387FF', '#55C667FF','#95D840FF', '#FDE725FF', 'gray50'), 
       lty = c(1,2,3,4,5,6,1,2,3,4), lwd = c(1,2,3,1,2,3,1,2,3,1), cex = 0.8, bty = "n")
dev.off()

#borderline vs the rest

#MLR
testdf$borderline <- as.factor(ifelse(testdf$outcome5 == 2, 1, 0))
m1borderline <- NB.imp(pred=MM1_pred_bo, outcome=borderline, data=testdf, imp=imp)$Results
m1borderline

#ridge 
m4borderline <- NB.imp(pred=MM4_pred_bo, outcome=borderline, data=testdf, imp=imp)$Results
m4borderline
#firth
m5borderline <- NB.imp(pred=MM5_pred_bo, outcome=borderline, data=testdf, imp=imp)$Results
m5borderline
#linear
m2borderline <- NB.imp(pred=MM2_pred_bo, outcome=borderline, data=testdf, imp=imp)$Results
m2borderline
#rf
m7borderline <- NB.imp(pred=MM7_pred_bo, outcome=borderline, data=testdf, imp=imp)$Results
m7borderline
#xgb
m8borderline <- NB.imp(pred=MM8_pred_bo, outcome=borderline, data=testdf, imp=imp)$Results
m8borderline
#nn
m9borderline <- NB.imp(pred=MM9_pred_bo, outcome=borderline, data=testdf, imp=imp)$Results
m9borderline
#svm
m10borderline <- NB.imp(pred=MM10_pred_bo, outcome=borderline, data=testdf, imp=imp)$Results
m10borderline


borderlinenb <- data.frame(m1borderline$CutOff,m1borderline$NB.mean, m4borderline$NB.mean, m5borderline$NB.mean, m2borderline$NB.mean, m7borderline$NB.mean, m8borderline$NB.mean, m9borderline$NB.mean, m10borderline$NB.mean, m1borderline$NBTA.mean)

#plot
tiff("netbenefit borderline.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(borderlinenb$m1borderline.CutOff,borderlinenb$m1borderline.NBTA.mean, type = "l", col = "#440154FF", lty=1, xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.5, 0.10), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(borderlinenb$m1borderline.CutOff, borderlinenb$m1borderline.NB.mean, type='l', col = "#482677FF", lty = 2, lwd = 2)
points(borderlinenb$m1borderline.CutOff, borderlinenb$m4borderline.NB.mean, type = 'l', col = "#404788FF", lty = 3, lwd = 3)
points(borderlinenb$m1borderline.CutOff, borderlinenb$m5borderline.NB.mean, type = 'l', col = "#33638DFF", lty = 4, lwd = 1)
points(borderlinenb$m1borderline.CutOff, borderlinenb$m2borderline.NB.mean, type = 'l', col = "#238A8DFF", lty=5,lwd = 2)
points(borderlinenb$m1borderline.CutOff, borderlinenb$m7borderline.NB.mean, type = 'l', col = "#20A387FF", lty=6,lwd = 3)
points(borderlinenb$m1borderline.CutOff, borderlinenb$m8borderline.NB.mean, type = 'l', col = "#55C667FF", lty = 1, lwd = 1)
points(borderlinenb$m1borderline.CutOff, borderlinenb$m9borderline.NB.mean, type = 'l', col = "#95D840FF", lty = 2, lwd = 2)
points(borderlinenb$m1borderline.CutOff, borderlinenb$m10borderline.NB.mean, type = 'l', col = "#FDE725FF", lty=4,lwd = 1)
abline(a = 0, b = 0, lty = 5, lwd = 1, col = "gray50")
legend("bottomleft", legend = c( "Treat all", "MLR", "Ridge MLR", "Firth MLR", "Linear MLR", "Random forest", "XGBoost", "Neural network", "Support vector machine", "Treat none"), ncol = 1, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("#440154FF", "#482677FF", "#404788FF",  "#33638DFF", "#238A8DFF", '#20A387FF', '#55C667FF','#95D840FF', '#FDE725FF', 'gray50'), 
       lty = c(1,2,3,4,5,6,1,2,3,4), lwd = c(1,2,3,1,2,3,1,2,3,1), cex = 0.8, bty = "n")
dev.off()

#Stage I vs rest
testdf$stage1 <- as.factor(ifelse(testdf$outcome5 == 3, 1, 0))

#MLR
m1stage1 <- NB.imp(pred=MM1_pred_bs1, outcome=stage1, data=testdf, imp=imp)$Results
m1stage1
#ridge
m4stage1 <- NB.imp(pred=MM4_pred_bs1, outcome=stage1, data=testdf, imp=imp)$Results
m4stage1
#firth
m5stage1 <- NB.imp(pred=MM5_pred_bs1, outcome=stage1, data=testdf, imp=imp)$Results
m5stage1
#lineqr
m2stage1 <- NB.imp(pred=MM2_pred_bs1, outcome=stage1, data=testdf, imp=imp)$Results
m2stage1
#rf
m7stage1 <- NB.imp(pred=MM7_pred_bs1, outcome=stage1, data=testdf, imp=imp)$Results
m7stage1
#xg
m8stage1 <- NB.imp(pred=MM8_pred_bs1, outcome=stage1, data=testdf, imp=imp)$Results
m8stage1
#nn
m9stage1 <- NB.imp(pred=MM9_pred_bs1, outcome=stage1, data=testdf, imp=imp)$Results
m9stage1
#svm
m10stage1 <- NB.imp(pred=MM10_pred_bs1, outcome=stage1, data=testdf, imp=imp)$Results
m10stage1


stage1nb <- data.frame(m1stage1$CutOff,m1stage1$NB.mean, m4stage1$NB.mean, m5stage1$NB.mean, m2stage1$NB.mean, m7stage1$NB.mean, m8stage1$NB.mean, m9stage1$NB.mean, m10stage1$NB.mean, m1stage1$NBTA.mean)
colnames(stage1nb) <- c("Cutoff", "m1", "m4", "m5","m2", "m7", "m8", "m9", "m10", "Treatall")

#plot
tiff("netbenefit stage1.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(stage1nb$Cutoff,stage1nb$Treatall, type = "l", col = "#440154FF", lty=1, xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.5, 0.10), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(stage1nb$Cutoff, stage1nb$m1, type='l', col = "#482677FF", lty = 2, lwd = 2)
points(stage1nb$Cutoff, stage1nb$m4, type = 'l', col = "#404788FF", lty = 3, lwd = 3)
points(stage1nb$Cutoff, stage1nb$m5, type = 'l', col = "#33638DFF", lty = 4, lwd = 1)
points(stage1nb$Cutoff, stage1nb$m2, type = 'l', col = "#238A8DFF", lty=5,lwd = 2)
points(stage1nb$Cutoff, stage1nb$m7, type = 'l', col = "#20A387FF", lty=6,lwd = 3)
points(stage1nb$Cutoff, stage1nb$m8, type = 'l', col = "#55C667FF", lty = 1, lwd = 1)
points(stage1nb$Cutoff, stage1nb$m9, type = 'l', col = "#95D840FF", lty = 2, lwd = 2)
points(stage1nb$Cutoff, stage1nb$m10, type = 'l', col = "#FDE725FF", lty=4,lwd = 1)
abline(a = 0, b = 0, lty = 5, lwd = 1, col = "gray50")
legend("bottomleft", legend = c( "Treat all", "MLR", "Ridge MLR", "Firth MLR", "Linear MLR", "Random forest", "XGBoost", "Neural network", "Support vector machine", "Treat none"), ncol = 1, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("#440154FF", "#482677FF", "#404788FF",  "#33638DFF", "#238A8DFF", '#20A387FF', '#55C667FF','#95D840FF', '#FDE725FF', 'gray50'), 
       lty = c(1,2,3,4,5,6,1,2,3,4), lwd = c(1,2,3,1,2,3,1,2,3,1), cex = 0.8, bty = "n")
dev.off()


#stage II-IV vs rest
testdf$stage24 <- as.factor(ifelse(testdf$outcome5 == 4, 1, 0))
m1stage24 <- NB.imp(pred=MM1_pred_bs2, outcome=stage24, data=testdf, imp=imp)$Results
m1stage24

#ridge
m4stage24 <- NB.imp(pred=MM4_pred_bs2, outcome=stage24, data=testdf, imp=imp)$Results
m4stage24
#firth
m5stage24 <- NB.imp(pred=MM5_pred_bs2, outcome=stage24, data=testdf, imp=imp)$Results
m5stage24
#lineqr
m2stage24 <- NB.imp(pred=MM2_pred_bs2, outcome=stage24, data=testdf, imp=imp)$Results
m2stage24
#rf
m7stage24 <- NB.imp(pred=MM7_pred_bs2, outcome=stage24, data=testdf, imp=imp)$Results
m7stage24
#xg
m8stage24 <- NB.imp(pred=MM8_pred_bs2, outcome=stage24, data=testdf, imp=imp)$Results
m8stage24
#nn
m9stage24 <- NB.imp(pred=MM9_pred_bs2, outcome=stage24, data=testdf, imp=imp)$Results
m9stage24
#svm
m10stage24 <- NB.imp(pred=MM10_pred_bs2, outcome=stage24, data=testdf, imp=imp)$Results
m10stage24


stage24nb <- data.frame(m1stage24$CutOff,m1stage24$NB.mean, m4stage24$NB.mean, m5stage24$NB.mean, m2stage24$NB.mean, m7stage24$NB.mean, m8stage24$NB.mean, m9stage24$NB.mean, m10stage24$NB.mean, m1stage24$NBTA.mean)
colnames(stage24nb) <- c("Cutoff", "m1", "m4", "m5","m2", "m7", "m8", "m9", "m10", "Treatall")

#plot
tiff("netbenefit stage24.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(stage24nb$Cutoff,stage24nb$Treatall, type = "l", col = "#440154FF", lty=1, xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.5, 0.20), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(stage24nb$Cutoff, stage24nb$m1, type='l', col = "#482677FF", lty = 2, lwd = 2)
points(stage24nb$Cutoff, stage24nb$m4, type = 'l', col = "#404788FF", lty = 3, lwd = 3)
points(stage24nb$Cutoff, stage24nb$m5, type = 'l', col = "#33638DFF", lty = 4, lwd = 1)
points(stage24nb$Cutoff, stage24nb$m2, type = 'l', col = "#238A8DFF", lty=5,lwd = 2)
points(stage24nb$Cutoff, stage24nb$m7, type = 'l', col = "#20A387FF", lty=6,lwd = 3)
points(stage24nb$Cutoff, stage24nb$m8, type = 'l', col = "#55C667FF", lty = 1, lwd = 1)
points(stage24nb$Cutoff, stage24nb$m9, type = 'l', col = "#95D840FF", lty = 2, lwd = 2)
points(stage24nb$Cutoff, stage24nb$m10, type = 'l', col = "#FDE725FF", lty=4,lwd = 1)
abline(a = 0, b = 0, lty = 5, lwd = 1, col = "gray50")
legend("bottomleft", legend = c( "Treat all", "MLR", "Ridge MLR", "Firth MLR", "Linear MLR", "Random forest", "XGBoost", "Neural network", "Support vector machine", "Treat none"), ncol = 1, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("#440154FF", "#482677FF", "#404788FF",  "#33638DFF", "#238A8DFF", '#20A387FF', '#55C667FF','#95D840FF', '#FDE725FF', 'gray50'), 
       lty = c(1,2,3,4,5,6,1,2,3,4), lwd = c(1,2,3,1,2,3,1,2,3,1), cex = 0.8, bty = "n")
dev.off()

#metastatic
testdf$meta <- as.factor(ifelse(testdf$outcome5 == 5, 1, 0))
m1meta <- NB.imp(pred=MM1_pred_mt, outcome=meta, data=testdf, imp=imp)$Results
m1meta

#ridge
m4meta <- NB.imp(pred=MM4_pred_mt, outcome=meta, data=testdf, imp=imp)$Results

#firth
m5meta <- NB.imp(pred=MM5_pred_mt, outcome=meta, data=testdf, imp=imp)$Results

#lineqr
m2meta <- NB.imp(pred=MM2_pred_mt, outcome=meta, data=testdf, imp=imp)$Results

#rf
m7meta <- NB.imp(pred=MM7_pred_mt, outcome=meta, data=testdf, imp=imp)$Results

#xg
m8meta <- NB.imp(pred=MM8_pred_mt, outcome=meta, data=testdf, imp=imp)$Results

#nn
m9meta <- NB.imp(pred=MM9_pred_mt, outcome=meta, data=testdf, imp=imp)$Results

#svm
m10meta <- NB.imp(pred=MM10_pred_mt, outcome=meta, data=testdf, imp=imp)$Results



metanb <- data.frame(m1meta$CutOff,m1meta$NB.mean, m4meta$NB.mean, m5meta$NB.mean, m2meta$NB.mean, m7meta$NB.mean, m8meta$NB.mean, m9meta$NB.mean, m10meta$NB.mean, m1meta$NBTA.mean)
colnames(metanb) <- c("Cutoff", "m1", "m4", "m5","m2", "m7", "m8", "m9", "m10", "Treatall")

#plot
tiff("netbenefit meta.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(metanb$Cutoff,metanb$Treatall, type = "l", col = "#440154FF", lty=1, xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.5, 0.10), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(metanb$Cutoff, metanb$m1, type='l', col = "#482677FF", lty = 2, lwd = 2)
points(metanb$Cutoff, metanb$m4, type = 'l', col = "#404788FF", lty = 3, lwd = 3)
points(metanb$Cutoff, metanb$m5, type = 'l', col = "#33638DFF", lty = 4, lwd = 1)
points(metanb$Cutoff, metanb$m2, type = 'l', col = "#238A8DFF", lty=5,lwd = 2)
points(metanb$Cutoff, metanb$m7, type = 'l', col = "#20A387FF", lty=6,lwd = 3)
points(metanb$Cutoff, metanb$m8, type = 'l', col = "#55C667FF", lty = 1, lwd = 1)
points(metanb$Cutoff, metanb$m9, type = 'l', col = "#95D840FF", lty = 2, lwd = 2)
points(metanb$Cutoff, metanb$m10, type = 'l', col = "#FDE725FF", lty=4,lwd = 1)
abline(a = 0, b = 0, lty = 5, lwd = 1, col = "gray50")
legend("bottomleft", legend = c( "Treat all", "MLR", "Ridge MLR", "Firth MLR", "Linear MLR", "Random forest", "XGBoost", "Neural network", "Support vector machine", "Treat none"), ncol = 1, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("#440154FF", "#482677FF", "#404788FF",  "#33638DFF", "#238A8DFF", '#20A387FF', '#55C667FF','#95D840FF', '#FDE725FF', 'gray50'), 
       lty = c(1,2,3,4,5,6,1,2,3,4), lwd = c(1,2,3,1,2,3,1,2,3,1), cex = 0.8, bty = "n")
dev.off()

