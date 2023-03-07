# Risk models for ovarian cancer diagnosis using different algorithms 

library(rms)
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
library(metamisc)
library(irr)
library(tidyverse)

            ## ------------- ##
            ## Without CA125 ##
            ## ------------- ##

# MLR

mm1 <- multinom(outcome5~oncocenter+Age+agercs+propsol+propsolrcs+papnrli+loc10+Shadows+Ascites+lesdmax+lesdmaxrcs, data=train)

# RIDGE MLR
#prepare data
y_train_std <- data.matrix(train[,"outcome5"])
x_train_std <- data.matrix(train[, c("lesdmax","lesdmaxrcs", "Age","agercs","oncocenter" ,"propsol","propsolrcs","papnrli", "loc10", "Shadows", "Ascites")])
y_test_std <- data.matrix(test[,"outcome5"])
x_test_std <- data.matrix(test[, c("lesdmax","lesdmaxrcs", "Age","agercs","oncocenter" ,"propsol","propsolrcs","papnrli", "loc10", "Shadows", "Ascites")])

control <- trainControl(method="cv", number=10, classProbs=TRUE, summaryFunction=multiClassSummary, search="grid")
lambda.param = 2^runif(30, min = -10, 3)
set.seed(234)
ridge<-train(y = make.names(y_train_std),
             x = x_train_std,
             method = 'glmnet',
             metric="logLoss",
             trControl=control,
             tuneGrid = expand.grid(alpha = 0, lambda = lambda.param)) 
ridge
mm4_pred <- ridge %>% predict(x_test_std, type="prob")

# RANDOM FOREST 
control <- trainControl(method="cv", number=10, classProbs=TRUE, summaryFunction=multiClassSummary, search="random")
set.seed(1234) 
mm7 <- train(make.names(outcome5)~oncocenter+Age+propsol+papnr+loc10+Shadows+Ascites+lesdmax, 
             data = train, 
             method = "ranger",
             trControl = control,
             tuneLength=30,
             metric="logLoss")
mm7

# XGBoost 
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


# NEURAL NETWORKS 
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


# SVM 
set.seed(789)
mm10 <- train(form=make.names(V1)~.,
                          method = "svmRadial",# Radial basis function kernel
                          tuneLength = 30,
                          metric = "logLoss",
                          data=d_training,
                          trControl=ctrl)  
mm10

###################################
#Boxplots of estimated probability#
###################################

mlr1 <- ggplot(mm1_pred1, aes(x=outcome55, y=X1))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), #-0.4
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))


mlr2 <- ggplot(mm1_pred1, aes(x=outcome55, y=X2))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))

mlr3 <- ggplot(mm1_pred1, aes(x=outcome55, y=X3))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(mm1_pred1, aes(x=outcome55, y=X4))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))

mlr5 <- ggplot(mm1_pred1, aes(x=outcome55, y=X5))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/MLR.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(), #+theme(axis.title.y=element_blank())
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')


grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob))

dev.off()



#ridge mlr

mlr1 <- ggplot(mm1_pred1, aes(x=outcome55, y=X1))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr2 <- ggplot(mm1_pred1, aes(x=outcome55, y=X2))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))

mlr3 <- ggplot(mm1_pred1, aes(x=outcome55, y=X3))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(mm1_pred1, aes(x=outcome55, y=X4))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))

mlr5 <- ggplot(mm1_pred1, aes(x=outcome55, y=X5))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/RIDGE MLR.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(),
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')
grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob))
dev.off()


#rf

mlr1 <- ggplot(mm1_pred1, aes(x=outcome55, y=X1))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))



mlr2 <- ggplot(mm1_pred1, aes(x=outcome55, y=X2))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(mm1_pred1, aes(x=outcome55, y=X3))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(mm1_pred1, aes(x=outcome55, y=X4))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr5 <- ggplot(mm1_pred1, aes(x=outcome55, y=X5))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/RF.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(),
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')

grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob))
dev.off()

#xg
mlr1 <- ggplot(mm1_pred1, aes(x=outcome55, y=X1))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))


mlr2 <- ggplot(mm1_pred1, aes(x=outcome55, y=X2))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(mm1_pred1, aes(x=outcome55, y=X3))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(mm1_pred1, aes(x=outcome55, y=X4))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr5 <- ggplot(mm1_pred1, aes(x=outcome55, y=X5))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/XGBoost.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(),
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')

grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob))
dev.off()

#nn

mlr1 <- ggplot(mm1_pred1, aes(x=outcome55, y=X1))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))
mlr1

mlr2 <- ggplot(mm1_pred1, aes(x=outcome55, y=X2))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(mm1_pred1, aes(x=outcome55, y=X3))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(mm1_pred1, aes(x=outcome55, y=X4))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr5 <- ggplot(mm1_pred1, aes(x=outcome55, y=X5))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/NN.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(),
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')

grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob))
dev.off()

#svm
mlr1 <- ggplot(mm1_pred1, aes(x=outcome55, y=X1))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))

mlr2 <- ggplot(mm1_pred1, aes(x=outcome55, y=X2))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(mm1_pred1, aes(x=outcome55, y=X3))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(mm1_pred1, aes(x=outcome55, y=X4))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr5 <- ggplot(mm1_pred1, aes(x=outcome55, y=X5))+
  geom_jitter(color="grey", size=0.2, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/SVM.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(),
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')

y.grob <- textGrob("Predicted probability of being ...", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)

x.grob <- textGrob("True outcome category", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))


grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob))

dev.off()

## ---------------------------------------- ##
##Scatterplots matrix of the estimated risks##
## -----------------------------------------##

#scatterplot of estimated risk of a benign tumour 
tiff("Descriptive/Scatterplot risk of benign.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(benigncombined[,1:6], pch=20, upper.panel=NULL, cex=0.9,  cex.labels=3, xlim=c(0,1), ylim=c(0,1)) 
dev.off()


#scatterplot of estimated risk of a borderline tumour 

tiff("Descriptive/Scatterplot borderline.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(bocombined[,1:6], pch=20, upper.panel=NULL, cex=0.9,  cex.labels=3, xlim=c(0,1), ylim=c(0,1)) 
dev.off()

#Scatterplot of estimated risk of a stage I primary invasive tumour

tiff("Descriptive/Scatterplot stage 1.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(s1combined[,1:6], pch=20, upper.panel=NULL, cex=0.9,  cex.labels=3, xlim=c(0,1), ylim=c(0,1))
dev.off()

#Scatterplot of estimated risk of a stage II-IV primary invasive tumour

tiff("Descriptive/Scatterplot stage 2 4.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(s2combined[,1:6], pch=20, upper.panel=NULL, cex=0.9,  cex.labels=3, xlim=c(0,1), ylim=c(0,1)) 
dev.off()

#Scatterplots of estimated risk of a secondary metastatic tumour

tiff("Descriptive/Scatterplot s meta.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(mtcombined[,1:6], pch=20, upper.panel=NULL, cex=0.9,  cex.labels=3, xlim=c(0,1), ylim=c(0,1))
dev.off()

## ----- ##
## Range ##
## ----- ##

range_benign <- apply(X=test[, c(172,177,182,187,192,197)], MARGIN=1, FUN=range)
test$Benign <- range_benign[2,]-range_benign[1,]
range_bord <- apply(X=test[, c(173,178,183,188,193,198)], MARGIN=1, FUN=range)
test$Borderline <- range_bord[2,]-range_bord[1,]
range_s1 <- apply(X=test[, c(174,179,184,189,194,199)], MARGIN=1, FUN=range)
test$'Stage I' <- range_s1[2,]-range_s1[1,]
range_s2 <- apply(X=test[, c(175,180,185,190,195,200)], MARGIN=1, FUN=range)
test$'Stage II-IV' <- range_s2[2,]-range_s2[1,]
range_smeta <- apply(X=test[, c(176,181,186,191,196,201)], MARGIN=1, FUN=range)
test$'Sec. metastatic' <- range_smeta[2,]-range_smeta[1,]

gg.woca <- dfwoca %>%
  ggplot(aes(x=Outcome, y=value))+
  geom_jitter(color="grey", size=0.1, alpha=0.8)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+ 
  #theme_ipsum()+
  xlab("Estimated probability")+
  ylab("Range of estimated probabilities across the models")+
  theme_light()+
  ylim(0,1)

test$MLR <- ifelse(test$mm1_pred>=0.10, "Malignant", "Benign")
test$Ridge <- ifelse(test$mm4_pred>=0.10, "Malignant", "Benign")
test$RF <- ifelse(test$mm7_pred>=0.10, "Malignant", "Benign")
test$XGBoost <- ifelse(test$mm8_pred>=0.10, "Malignant", "Benign")
test$NN <- ifelse(test$mm9_pred>=0.10, "Malignant", "Benign")
test$SVM <- ifelse(test$mm10_pred>=0.10, "Malignant", "Benign")

ratings <- test %>% select(Ridge,SVM)
100-as.numeric(agree(ratings)[5])



#Combine the smallest centres into one
levels(test$Center)[levels(test$Center)=="IUK"] <-"Other"
levels(test$Center)[levels(test$Center)=="MIT"] <-"Other"
levels(test$Center)[levels(test$Center)=="FLI"] <-"Other"
levels(test$Center)[levels(test$Center)=="NUK"] <-"Other"
levels(test$Center)[levels(test$Center)=="FIT"] <-"Other"
levels(test$Center)[levels(test$Center)=="CEG"] <-"Other"

#############################
##1. AUC BENIGN VS MALIGNANT#
#############################
#MLR
mm1_extauc <- AUC.IOTA.BAYES(pred=mm1_pred, outcome=outcome1, center=Center, data=test)

#Ridge
mm2_extauc <- AUC.IOTA.BAYES(pred=mm4_pred, outcome=outcome1, center=Center, data=test)

#RF
mm5_extauc <- AUC.IOTA.BAYES(pred=mm7_pred, outcome=outcome1, center=Center, data=test)

#XGBOOST
mm6_extauc <- AUC.IOTA.BAYES(pred=mm8_pred, outcome=outcome1, center=Center, data=test)

##NN
mm7_extauc <- AUC.IOTA.BAYES(pred=mm9_pred, outcome=outcome1, center=Center, data=test)

#SVM
mm8_extauc <- AUC.IOTA.BAYES(pred=mm10_pred, outcome=outcome1, center=Center, data=test)

#Overall forest plot
NA.forest <- mm1_extauc$Performance[1,]
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, mm1_extauc$Performance[1,],mm2_extauc$Performance[1,],
                     mm5_extauc$Performance[1,], mm6_extauc$Performance[1,], mm7_extauc$Performance[1,], mm8_extauc$Performance[1,])
Summary.AUC$Model <- c('', 'Standard MLR', 'Ridge MLR',  'Random forest', 'XGBoost','Neural network', 'Support vector machine')

Summary.PI <- rbind(NA.forest, mm1_extauc$Performance[2,], mm2_extauc$Performance[2,],
                    mm5_extauc$Performance[2,], mm6_extauc$Performance[2,],mm7_extauc$Performance[2,],mm8_extauc$Performance[2,])
Summary.PI$Model <- c('', 'Standard MLR', 'Ridge MLR',  'Random forest', 'XGBoost','Neural network', 'Support vector machine')

Summary.AUC
tabletext <- cbind(
  c('Model', 'Standard MLR', 'Ridge MLR', 'Random forest', 'XGBoost','Neural network', 'Support vector machine'),
  c('AUROC (95% CI)', 
    paste(format(round(mm1_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm1_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm1_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm2_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm2_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm2_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm5_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm5_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm5_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm6_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm6_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm6_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm7_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm7_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm7_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm8_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm8_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm8_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(mm1_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm1_extauc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm2_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm2_extauc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm5_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm5_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm6_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm6_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm7_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm7_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm8_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm8_extauc$Performance$UL[2], 2), nsmall = 2), ")"))
)


tiff("Binary AUC/Summary plot.tiff", width = 31, height = 20, units = "cm", res = 300) #20 groups rcs/

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
           xticks = c(0.80, 0.9, 1), xlog = TRUE, clip = c(0.8, 1),
           zero=NA)
dev.off()

#### CLINICAL UTILITY ####

#MLR
CULR1<- DataWinBugs(pred=test$mm1_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`[[10]]

#linear MLR
CULR2<- DataWinBugs(pred=test$mm2_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`[[10]]

#ridge
CULR4<- DataWinBugs(pred=test$mm4_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`[[10]]

#firth
CULR5<- DataWinBugs(pred=test$mm5_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`[[10]]

#rf
CULR7<- DataWinBugs(pred=test$mm7_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`[[10]]

#xg
CULR8 <- DataWinBugs(pred=test$mm8_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`[[10]]

#nn
CULR9<- DataWinBugs(pred=test$mm9_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`[[10]]

#svm
CULR10<- DataWinBugs(pred=test$mm10_pred, outcome=test$outcome1, center=test$Center, data=test)$`Results`[[10]]

# plot

tiff("NB/nb without ca125.tiff", width = 15, height = 15, units = "cm", res = 300)
plot(datanb$threshold,datanb$`NB treatall`, type = "l", col = "gray50", lty=1, xlab = "Threshold probability of any malignancy", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.01, 0.35), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(datanb$threshold, datanb$`NB LR1`, type='l', col = "#3399CC", lty = 2, lwd = 2)
points(datanb$threshold,datanb$`NB RIDGE`, type = 'l', col = "#CCFF66", lty = 3, lwd = 3)
points(datanb$threshold, datanb$`NB RF`, type = 'l', col = "#332288", lty = 6, lwd = 3)
points(datanb$threshold, datanb$`NB XGB`, type = 'l', col = "#00FF66", lty=1,lwd = 1)
points(datanb$threshold, datanb$`NB NN`, type = 'l', col = "#FF33CC", lty=2,lwd = 2)
points(datanb$threshold, datanb$`NB SVM`, type='l', col = "#882255", lty = 3, lwd = 3)
abline(a = 0, b = 0, lty = 4, lwd = 1, col = "gray50")
legend(x=0.05, y=0.05,legend = c("Standard MLR", "Ridge MLR",   "Random forest", "XGBoost", "Neural network", "Support vector machine","Treat all","Treat none"), ncol = 3, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c( "#3399CC", "#CCFF66",   '#332288', '#00FF66','#FF33CC',"#882255","gray50", 'gray50'), 
       lty = c(2,3,6,1,2,3,1,4), lwd = c(2,3,3,1,2,3,1,1), cex = 0.68, bty = "n")
dev.off()

######
#Combine together smallest centers to calculate PDI and calibration#

levels(test$Center)[levels(test$Center)=="IUK"] <-"Other"
levels(test$Center)[levels(test$Center)=="MIT"] <-"Other"
levels(test$Center)[levels(test$Center)=="FLI"] <-"Other"
levels(test$Center)[levels(test$Center)=="NUK"] <-"Other"
levels(test$Center)[levels(test$Center)=="FIT"] <-"Other"
levels(test$Center)[levels(test$Center)=="CEG"] <-"Other"
levels(test$Center)[levels(test$Center)=="SIT"] <-"Other"
levels(test$Center)[levels(test$Center)=="UDI"] <-"Other"
levels(test$Center)[levels(test$Center)=="BAI"] <-"Other"
levels(test$Center)[levels(test$Center)=="LPO"] <-"Other"
levels(test$Center)[levels(test$Center)=="NCI"] <-"Other"
levels(test$Center)[levels(test$Center)=="PSP"] <-"Other"
levels(test$Center)[levels(test$Center)=="TIT"] <-"Other"
levels(test$Center)[levels(test$Center)=="MPO"] <-"Other"
levels(test$Center)[levels(test$Center)=="BCH"] <-"Other"

#########
## PDI ##
#########

#MLR
m1extpdi <- PDI.centre(pred=mm1_pred, outcome=outcome5, center=Center,  data=test)

#Ridge
m2extpdi <- PDI.centre(pred=mm4_pred, outcome=outcome5, center=Center,  data=test)

#RF
m5extpdi <- PDI.centre(pred=mm7_pred, outcome=outcome5, center=Center, data=test)

#XGBoost
m6extpdi <- PDI.centre(pred=mm8_pred, outcome=outcome5, center=Center, data=test)

#NN
m7extpdi <- PDI.centre(pred=mm9_pred, outcome=outcome5, center=Center,  data=test)

#SVM
m8extpdi <- PDI.centre(pred=mm10_pred, outcome=outcome5, center=Center, data=test)


#Overall forest plot
NA.forest <- m1extpdi$Performance[1,]
NA.forest <- NA

Summary.PDI <- rbind(NA.forest, m1extpdi$Performance[1,],m2extpdi$Performance[1,],  m5extpdi$Performance[1,],
                     m6extpdi$Performance[1,],m7extpdi$Performance[1,],m8extpdi$Performance[1,])
Summary.PDI$Model <- c('',  'Standard MLR', 'Ridge MLR',  "Random forest", "XGBoost", "Neural network", "Support vector machine")

Summary.PI <- rbind(NA.forest, m1extpdi$Performance[2,], m2extpdi$Performance[2,], 
                    m5extpdi$Performance[2,],m6extpdi$Performance[2,],
                    m7extpdi$Performance[2,],m8extpdi$Performance[2,]
)
Summary.PI$Model <- c('',  'Standard MLR', 'Ridge MLR',  "Random forest", "XGBoost", "Neural network", "Support vector machine")


tabletext <- cbind(
  c('Model',  'Standard MLR', 'Ridge MLR',  "Random forest", "XGBoost", "Neural network", "Support vector machine"),
  c('PDI (95% CI)', 
    paste(format(round(m1extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m1extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m1extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m2extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m2extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m2extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m5extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m5extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m5extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m6extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m6extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m6extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m7extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m7extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m7extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m8extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m8extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m8extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = "")
  ),
  c('95% PI', 
    paste0("(", format(round(m1extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m1extpdi$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(m2extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m2extpdi$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(m5extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m5extpdi$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(m6extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m6extpdi$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(m7extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m7extpdi$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(m8extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m8extpdi$Performance$UL[2], 2), nsmall = 2), ")")
  )
)


tiff("PDI/summary plot.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.PDI$PDI, 3),
           lower = round(Summary.PDI$LL, 3),
           upper = round(Summary.PDI$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE, TRUE),
           xlab = "PDI (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           zero=NA,
           xticks = c(0.40,0.50,0.60, 0.70), xlog = TRUE, clip = c(0.40, 70))
dev.off()

######################
#### PAIRWISE AUC ####
######################

##Pair 1 benign vs borderline##
#MLR
test$bvsbridge <- mm1_pred[, 1] / (mm1_pred[,1] + mm1_pred[,2])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==2,])

#Ridge
test$bvsbridge <- mm4_pred[, 1] / (mm4_pred[,1] + mm4_pred[,2])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==2,])

#RF
test$bvsbridge <- mm7_pred[, 1] / (mm7_pred[,1] + mm7_pred[,2])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==2,])

#XGBoost
test$bvsbridge <- mm8_pred[, 1] / (mm8_pred[,1] + mm8_pred[,2])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==2,])

#NN
test$bvsbridge <- mm9_pred[, 1] / (mm9_pred[,1] + mm9_pred[,2])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==2,])

#SVM
test$bvsbridge <- mm10_pred[, 1] / (mm10_pred[,1] + mm10_pred[,2])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==1|test$outcome5==2,])

##Pair 2 benign vs stage 1##
#MLR
test$bvsbridge <- mm1_pred[, 1] / (mm1_pred[,1] + mm1_pred[,3])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==1|test$outcome5==3,])

#Ridge
test$bvsbridge <- mm4_pred[, 1] / (mm4_pred[,1] + mm4_pred[,3])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==3,])

#RF
test$bvsbridge <- mm7_pred[, 1] / (mm7_pred[,1] + mm7_pred[,3])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==3,])

#XGBoost
test$bvsbridge <- mm8_pred[, 1] / (mm8_pred[,1] + mm8_pred[,3])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==1|test$outcome5==3,])

#NN
test$bvsbridge <- mm9_pred[, 1] / (mm9_pred[,1] + mm9_pred[,3])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==1|test$outcome5==3,])

#SVM
test$bvsbridge <- mm10_pred[, 1] / (mm10_pred[,1] + mm10_pred[,3])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==3,])


##Pair 3 benign vs stage 2-4##
#MLR
test$bvsbridge <- mm1_pred[, 1] / (mm1_pred[,1] + mm1_pred[,4])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==1|test$outcome5==4,])

#Ridge
test$bvsbridge <- mm4_pred[, 1] / (mm4_pred[,1] + mm4_pred[,4])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==4,])

#RF
test$bvsbridge <- mm7_pred[, 1] / (mm7_pred[,1] + mm7_pred[,4])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==4,])

#XGBoost
test$bvsbridge <- mm8_pred[, 1] / (mm8_pred[,1] + mm8_pred[,4])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==4,])

#NN
test$bvsbridge <- mm9_pred[, 1] / (mm9_pred[,1] + mm9_pred[,4])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==4,])

#SVM
test$bvsbridge <- mm10_pred[, 1] / (mm10_pred[,1] + mm10_pred[,4])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==1|test$outcome5==4,])

##Pair 4 benign vs secondary metastatic##
#MLR
test$bvsbridge <- mm1_pred[, 1] / (mm1_pred[,1] + mm1_pred[,5])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==5,])

#Ridge
test$bvsbridge <- mm4_pred[, 1] / (mm4_pred[,1] + mm4_pred[,5])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==5,])

#RF
test$bvsbridge <- mm7_pred[, 1] / (mm7_pred[,1] + mm7_pred[,5])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==1|test$outcome5==5,])

#XGBoost
test$bvsbridge <- mm8_pred[, 1] / (mm8_pred[,1] + mm8_pred[,5])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==5,])

#NN
test$bvsbridge <- mm9_pred[, 1] / (mm9_pred[,1] + mm9_pred[,5])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==1|test$outcome5==5,])

#SVM
test$bvsbridge <- mm10_pred[, 1] / (mm10_pred[,1] + mm10_pred[,5])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==1|test$outcome5==5,])


##Pair 5##
#MLR
test$bvsbridge <- mm1_pred[, 2] / (mm1_pred[,2] + mm1_pred[,3])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==2|test$outcome5==3,])

#Ridge
test$bvsbridge <- mm4_pred[, 2] / (mm4_pred[,2] + mm4_pred[,3])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==3,])

#RF
test$bvsbridge <- mm7_pred[, 2] / (mm7_pred[,2] + mm7_pred[,3])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==2|test$outcome5==3,])

#XGBoost
test$bvsbridge <- mm8_pred[, 2] / (mm8_pred[,2] + mm8_pred[,3])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==3,])

#NN
testdf$bvsbridge <- mm9_pred[, 2] / (mm9_pred[,2] + mm9_pred[,3])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==2|test$outcome5==3,])

#SVM
test$bvsbridge <- mm10_pred[, 2] / (mm10_pred[,2] + mm10_pred[,3])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==2|test$outcome5==3,])


##Pair 6##
#MLR
test$bvsbridge <- mm1_pred[, 2] / (mm1_pred[,2] + mm1_pred[,4])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==4,])

#Ridge
test$bvsbridge <- mm4_pred[, 2] / (mm4_pred[,2] + mm4_pred[,4])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==2|test$outcome5==4,])

#RF
test$bvsbridge <- mm7_pred[, 2] / (mm7_pred[,2] + mm7_pred[,4])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==4,])

#XGBoost
test$bvsbridge <- mm8_pred[, 2] / (mm8_pred[,2] + mm8_pred[,4])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==4,])

#NN
test$bvsbridge <- mm9_pred[, 2] / (mm9_pred[,2] + mm9_pred[,4])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==2|test$outcome5==4,])

#SVM
test$bvsbridge <- mm10_pred[, 2] / (mm10_pred[,2] + mm10_pred[,4])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==4,])


##Pair 7##
#MLR
test$bvsbridge <- mm1_pred[, 2] / (mm1_pred[,2] + mm1_pred[,5])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==2|test$outcome5==5,])

#Ridge
test$bvsbridge <- mm4_pred[, 2] / (mm4_pred[,2] + mm4_pred[,5])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==5,])

#RF
test$bvsbridge <- mm7_pred[, 2] / (mm7_pred[,2] + mm7_pred[,5])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==5,])

#XGBoost
test$bvsbridge <- mm8_pred[, 2] / (mm8_pred[,2] + mm8_pred[,5])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==5,])

#NN
test$bvsbridge <- mm9_pred[, 2] / (mm9_pred[,2] + mm9_pred[,5])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==5,])

#SVM
test$bvsbridge <- mm10_pred[, 2] / (mm10_pred[,2] + mm10_pred[,5])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==2|test$outcome5==5,])


##Pair 8##
#MLR
test$bvsbridge <- mm1_pred[, 3] / (mm1_pred[,3] + mm1_pred[,4])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==3|test$outcome5==4,])

#Ridge
test$bvsbridge <- mm4_pred[, 3] / (mm4_pred[,3] + mm4_pred[,4])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==4,])

#RF
test$bvsbridge <- mm7_pred[, 3] / (mm7_pred[,3] + mm7_pred[,4])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==4,])

#XGBoost
test$bvsbridge <- mm8_pred[, 3] / (mm8_pred[,3] + mm8_pred[,4])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==4,])

#NN
test$bvsbridge <- mm9_pred[, 3] / (mm9_pred[,3] + mm9_pred[,4])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==4,])

#SVM
test$bvsbridge <- mm10_pred[, 3] / (mm10_pred[,3] + mm10_pred[,4])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==4,])


##Pair 9##
#MLR
test$bvsbridge <- mm1_pred[, 3] / (mm1_pred[,3] + mm1_pred[,5])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==5,])

#Ridge
test$bvsbridge <- mm4_pred[, 3] / (mm4_pred[,3] + mm4_pred[,5])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==5,])

#RF
test$bvsbridge <- mm7_pred[, 3] / (mm7_pred[,3] + mm7_pred[,5])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==5,])

#XGBoost
test$bvsbridge <- mm8_pred[, 3] / (mm8_pred[,3] + mm8_pred[,5])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==5,])

#NN
test$bvsbridge <- mm9_pred[, 3] / (mm9_pred[,3] + mm9_pred[,5])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==3|test$outcome5==5,])

#SVM
test$bvsbridge <- mm10_pred[, 3] / (mm10_pred[,3] + mm10_pred[,5])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,data=test[test$outcome5==3|test$outcome5==5,])


##Pair 10##
#MLR
test$bvsbridge <- mm1_pred[, 4] / (mm1_pred[,4] + mm1_pred[,5])
bvsbAUC1 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==4|test$outcome5==5,])

#Ridge
test$bvsbridge <- mm4_pred[, 4] / (mm4_pred[,4] + mm4_pred[,5])
bvsbAUC2 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==4|test$outcome5==5,])

#RF
test$bvsbridge <- mm7_pred[, 4] / (mm7_pred[,4] + mm7_pred[,5])
bvsbAUC5 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==4|test$outcome5==5,])

#XGBoost
test$bvsbridge <- mm8_pred[, 4] / (mm8_pred[,4] + mm8_pred[,5])
bvsbAUC6 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==4|test$outcome5==5,])

#NN
test$bvsbridge <- mm9_pred[, 4] / (mm9_pred[,4] + mm9_pred[,5])
bvsbAUC7 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==4|test$outcome5==5,])

#SVM
test$bvsbridge <- mm10_pred[, 4] / (mm10_pred[,4] + mm10_pred[,5])
bvsbAUC8 <- AUC.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center, data=test[test$outcome5==4|test$outcome5==5,])


###-----------------------------------------------###
###--------------- With CA125 --------------------###
###-----------------------------------------------###

NRImp <- length(unique(testdf$imp))

#1. MLR#

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  tmp = traindf[traindf$imp == i, ]
  MM1 = multinom(outcome5~oncocenter+Age+agercs+propsol+propsolrcs+papnrli+loc10+Shadows+Ascites+lesdmax+lesdmaxrcs+logca125imp+logca125imprcs, data = tmp)
  
  
  return(MM1)
}, simplify = F, USE.NAMES = T)

ResultsMM1 = do.call("sapply", args = ArgzCalc)

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

tiff("MLR.tiff", width=14, height=14, units="cm", res=300)
mlr.hist <- ggplot(testdf, aes(x=MM1_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D"))+scale_y_continuous(breaks=c(-1000,0, 1000, 2000, 3000),labels = c(1000,0,1000,2000,3000))
mlr.hist
dev.off()


#Ridge logistic regression

control <- trainControl(method="cv", number=10, classProbs=TRUE, summaryFunction=multiClassSummary, search="grid")
set.seed(234)
lambda.param = 2^runif(30, min = -10, 3)

ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
  y_train_std <- data.matrix(traindf[traindf$imp == i,"outcome5"])
  x_train_std <- data.matrix(traindf[traindf$imp == i, c("lesdmax","lesdmaxrcs","agercs", "Age","oncocenter" ,"propsol","propsolrcs","papnrli", "loc10", "Shadows", "Ascites", "logca125imp", "logca125imprcs")])
  set.seed(234)
  MM4<-train(y = make.names(y_train_std),
             x = x_train_std,
             method = 'glmnet',
             metric="logLoss",
             trControl=control,
             tuneGrid = expand.grid(alpha = 0, lambda = lambda.param)) 
  return(MM4)
}, simplify = F, USE.NAMES = T)


ResultsMM4 = do.call("sapply", args = ArgzCalc)


#prediction on validation data

y_test_std <- data.matrix(testdf[,"outcome5"])
x_test_std <- data.matrix(testdf[, c("lesdmax","lesdmaxrcs", "Age","agercs","oncocenter" ,"propsol","propsolrcs","papnrli", "loc10", "Shadows", "Ascites", "logca125imp", "logca125imprcs")])

pMM4 <- data.frame(predict(ResultsMM4, x_test_std, type="prob"))


testdf$MM4_pred_be <- rowMeans(subset(pMM4, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
testdf$MM4_pred_bo <- rowMeans(subset(pMM4, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
testdf$MM4_pred_bs1 <- rowMeans(subset(pMM4, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
testdf$MM4_pred_bs2 <- rowMeans(subset(pMM4, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
testdf$MM4_pred_mt <- rowMeans(subset(pMM4, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)


testdf$MM4_mal <- rowSums(subset(testdf, select=c(MM4_pred_bo, MM4_pred_bs1, MM4_pred_bs2, MM4_pred_mt)))

MM4pred <- data.frame(testdf$MM4_pred_be, testdf$MM4_pred_bo, testdf$MM4_pred_bs1, testdf$MM4_pred_bs2, testdf$MM4_pred_mt)

tiff("Ridge MLR.tiff", width=14, height=14, units="cm", res=300)
ridgemlr.hist <- ggplot(testdf, aes(x=MM4_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(breaks = c(-1000,0,1000,2000, 3000), labels=c(1000,0,1000,2000,3000))
ridgemlr.hist
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

ResultsMM7 = do.call("sapply", args = ArgzCalc)

#prediction on validation data

pMM7 <- data.frame(predict(ResultsMM7, newdata=testdf, type="prob"))

testdf$MM7_pred_be <- rowMeans(subset(pMM7, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
testdf$MM7_pred_bo <- rowMeans(subset(pMM7, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
testdf$MM7_pred_bs1 <- rowMeans(subset(pMM7, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
testdf$MM7_pred_bs2 <- rowMeans(subset(pMM7, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9 )), na.rm = TRUE)
testdf$MM7_pred_mt <- rowMeans(subset(pMM7, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9 )), na.rm = TRUE)

testdf$MM7_mal <- rowSums(subset(testdf, select=c(MM7_pred_bo, MM7_pred_bs1, MM7_pred_bs2, MM7_pred_mt)))

MM7pred <- data.frame(testdf$MM7_pred_be, testdf$MM7_pred_bo, testdf$MM7_pred_bs1, testdf$MM7_pred_bs2, testdf$MM7_pred_mt)

tiff("RF.tiff", width=14, height=14, units="cm", res=300)
rf.hist <- ggplot(testdf, aes(x=MM7_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(breaks = c(-1000,0,1000,2000, 3000), labels=c(1000,0,1000,2000,3000))
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

tiff("XGBoost.tiff", width=14, height=14, units="cm", res=300)
xg.hist <- ggplot(testdf, aes(x=MM8_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(breaks = c(-1000,0,1000,2000, 3000, 4000, 5000), labels=c(1000,0,1000,2000,3000, 4000, 5000))
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

tiff("NN.tiff", width=14, height=14, units="cm", res=300)
nn.hist <- ggplot(testdf, aes(x=MM9_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(breaks = c(-1000,0,1000,2000, 3000, 4000, 5000), labels=c(1000,0,1000,2000,3000, 4000, 5000))
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


#prediction on test data

pMM10 <- data.frame(predict(ResultsMM10, newdata=d_test, type="prob"))

testdf$MM10_pred_be <- rowMeans(subset(pMM10, select = c(X1, X1.1,X1.2, X1.3,X1.4,X1.5,X1.6,X1.7,X1.8,X1.9 )), na.rm = TRUE)
testdf$MM10_pred_bo <- rowMeans(subset(pMM10, select = c(X2, X2.1,X2.2, X2.3,X2.4,X2.5,X2.6,X2.7,X2.8,X2.9 )), na.rm = TRUE)
testdf$MM10_pred_bs1 <- rowMeans(subset(pMM10, select = c(X3, X3.1,X3.2, X3.3,X3.4,X3.5,X3.6,X3.7,X3.8,X3.9 )), na.rm = TRUE)
testdf$MM10_pred_bs2 <- rowMeans(subset(pMM10, select = c(X4, X4.1,X4.2, X4.3,X4.4,X4.5,X4.6,X4.7,X4.8,X4.9)), na.rm = TRUE)
testdf$MM10_pred_mt <- rowMeans(subset(pMM10, select = c(X5, X5.1,X5.2, X5.3,X5.4,X5.5,X5.6,X5.7,X5.8,X5.9)), na.rm = TRUE)

testdf$MM10_mal <- rowSums(subset(testdf, select=c(MM10_pred_bo, MM10_pred_bs1, MM10_pred_bs2, MM10_pred_mt)))

MM10pred <- data.frame(testdf$MM10_pred_be, testdf$MM10_pred_bo, testdf$MM10_pred_bs1, testdf$MM10_pred_bs2, testdf$MM10_pred_mt)

tiff("SVM.tiff", width=14, height=14, units="cm", res=300)
svm.hist <- ggplot(testdf, aes(x=MM10_mal)) +
  geom_histogram(data=subset(testdf,outcome1 == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
  geom_histogram(data=subset(testdf,outcome1 == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
  theme_minimal() +
  labs(x = "Estimated risk of malignancy",
       y = "Frequency",
       fill = "Outcome") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
  scale_y_continuous(breaks = c(-1000,0,1000,2000, 3000, 4000, 5000), labels=c(1000,0,1000,2000,3000, 4000, 5000))
svm.hist

dev.off()

##############
#Scatterplots#
##############

#scatterplots of the probability of a benign tumour# 

tiff("Descriptive/Scatterplot risk of benign.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(benigncombined[,1:6], pch=".", upper.panel=NULL, cex=1.5,  cex.labels=3, xlim=c(0,1), ylim=c(0,1))
dev.off()


#scatterplot of estimated risk of a borderline tumour 

tiff("Descriptive/Scatterplot risk of borderline.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(bocombined[,1:6], pch=".", upper.panel=NULL, cex=1.5,  cex.labels=3, xlim=c(0,1), ylim=c(0,1))
dev.off()


#Scatterplot of estimated risk of a stage I primary invasive tumour

tiff("Descriptive/Scatterplot risk of stage 1.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(s1combined[,1:6], pch=".", upper.panel=NULL, cex=1.5,  cex.labels=3, xlim=c(0,1), ylim=c(0,1))
dev.off()

#Scatterplot of estimated risk of a stage II-IV primary invasive tumour

tiff("Descriptive/Scatterplot risk of stage 24.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(s2combined[,1:6], pch=".", upper.panel=NULL, cex=1.5,  cex.labels=3, xlim=c(0,1), ylim=c(0,1))
dev.off()

#Scatterplots of estimated risk of a secondary metastatic tumour

tiff("Descriptive/Scatterplot risk of smeta.tiff", width=30, height=30, units="cm", res=300)
par(cex.axis=1.2)
pairs(mtcombined[,1:6], pch=".", upper.panel=NULL, cex=1.5,  cex.labels=3, xlim=c(0,1), ylim=c(0,1))
dev.off()

#############################
#Boxplots of estimated risks#
#############################

mlr1 <- ggplot(testdf, aes(x=outcome55, y=MM1_pred_be))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), #-0.4
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))


mlr2 <- ggplot(testdf, aes(x=outcome55, y=MM1_pred_bo))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(testdf, aes(x=outcome55, y=MM1_pred_bs1))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(testdf, aes(x=outcome55, y=MM1_pred_bs2))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr5 <- ggplot(testdf, aes(x=outcome55, y=MM1_pred_mt))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/MLR.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(), #+theme(axis.title.y=element_blank())
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')


grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob)) #library(gridExtra)

dev.off()


#ridge mlr
mlr1 <- ggplot(testdf, aes(x=outcome55, y=MM4_pred_be))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), #-0.4
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))


mlr2 <- ggplot(testdf, aes(x=outcome55, y=MM4_pred_bo))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(testdf, aes(x=outcome55, y=MM4_pred_bs1))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(testdf, aes(x=outcome55, y=MM4_pred_bs2))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr5 <- ggplot(testdf, aes(x=outcome55, y=MM4_pred_mt))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/RIDGE.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(), #+theme(axis.title.y=element_blank())
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')

grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob)) 
dev.off()

#rf
mlr1 <- ggplot(testdf, aes(x=outcome55, y=MM7_pred_be))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), #-0.4
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))


mlr2 <- ggplot(testdf, aes(x=outcome55, y=MM7_pred_bo))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(testdf, aes(x=outcome55, y=MM7_pred_bs1))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(testdf, aes(x=outcome55, y=MM7_pred_bs2))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))

mlr5 <- ggplot(testdf, aes(x=outcome55, y=MM7_pred_mt))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/RF.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(), #+theme(axis.title.y=element_blank())
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')
y.grob <- textGrob("Predicted probability of being ...", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)

x.grob <- textGrob("True outcome category", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))

grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob)) 
dev.off()

#xg
mlr1 <- ggplot(testdf, aes(x=outcome55, y=MM8_pred_be))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), #-0.4
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))


mlr2 <- ggplot(testdf, aes(x=outcome55, y=MM8_pred_bo))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(testdf, aes(x=outcome55, y=MM8_pred_bs1))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))

mlr4 <- ggplot(testdf, aes(x=outcome55, y=MM8_pred_bs2))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr5 <- ggplot(testdf, aes(x=outcome55, y=MM8_pred_mt))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/XGBoost.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(), #+theme(axis.title.y=element_blank())
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')


grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob)) 
dev.off()

#nn
mlr1 <- ggplot(testdf, aes(x=outcome55, y=MM9_pred_be))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), #-0.4
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))


mlr2 <- ggplot(testdf, aes(x=outcome55, y=MM9_pred_bo))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(testdf, aes(x=outcome55, y=MM9_pred_bs1))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(testdf, aes(x=outcome55, y=MM9_pred_bs2))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr5 <- ggplot(testdf, aes(x=outcome55, y=MM9_pred_mt))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/NN.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(), #+theme(axis.title.y=element_blank())
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')


grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob)) 
dev.off()

#svm
mlr1 <- ggplot(testdf, aes(x=outcome55, y=MM10_pred_be))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+#xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.55, y=1.03, label="... Benign", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), #-0.4
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y=element_text(size=15))


mlr2 <- ggplot(testdf, aes(x=outcome55, y=MM10_pred_bo))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=1.6, y=1.03, label="... Borderline", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))


mlr3 <- ggplot(testdf, aes(x=outcome55, y=MM10_pred_bs1))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Stage I primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr4 <- ggplot(testdf, aes(x=outcome55, y=MM10_pred_bs2))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.9, y=1.03, label="... Stage II-IV primary invasive", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


mlr5 <- ggplot(testdf, aes(x=outcome55, y=MM10_pred_mt))+
  geom_jitter(color="grey", size=0.07, alpha=0.7)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+
  theme(legend.position = "none", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12))+
  #xlab("True class")+
  ylab("Predicted probability of being ...") + geom_text(x=2.8, y=1.03, label="... Secondary metastasis", size=4.5)+scale_y_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1))


tiff("Descriptive/SVM.tiff", width=31, height=26, units="cm", res=300)
plott <- cowplot::plot_grid(mlr1+theme(axis.title.y=element_blank()), mlr2+ theme(axis.text.y = element_blank(), #+theme(axis.title.y=element_blank())
                                                                                  axis.ticks.y = element_blank(),
                                                                                  axis.title.y = element_blank(),
                                                                                  axis.line.y = element_blank()), 
                            mlr3+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr4+ theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank(),
                                        axis.line.y = element_blank()), 
                            mlr5+ 
                              theme(axis.text.y = element_blank(),
                                    axis.line.y = element_blank(),
                                    axis.title.y= element_blank(),
                                    axis.ticks.y= element_blank()),
                            nrow = 1, 
                            align = 'h', axis = 'tb')

grid.arrange(arrangeGrob(plott, left = y.grob,bottom = x.grob)) 
dev.off()

############################
## AUC benign vs malignant##
############################
#Combine smallest centers
levels(testdf$Center)[levels(testdf$Center)=="IUK"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="MIT"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="FLI"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="NUK"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="FIT"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="CEG"] <-"Other"

### AUROC for benign vs malignant tumor
#MLR
mm1_extauc <- AUCimp.IOTA.BAYES(pred=MM1_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)

#RIDGE
mm4_extauc <- AUCimp.IOTA.BAYES(pred=MM4_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)

#RF
mm7_extauc <- AUCimp.IOTA.BAYES(pred=MM7_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)

#XG
mm8_extauc <- AUCimp.IOTA.BAYES(pred=MM8_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)

#NN
mm9_extauc <- AUCimp.IOTA.BAYES(pred=MM9_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)

#SVM
mm10_extauc <- AUCimp.IOTA.BAYES(pred=MM10_mal, outcome=outcome1, center=Center, imp=imp, data=testdf)

#Forest plot for every model

NA.forest <- mm1_extauc$Performance[1,]
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, mm1_extauc$Performance[1,],mm4_extauc$Performance[1,], 
                     mm7_extauc$Performance[1,], mm8_extauc$Performance[1,], mm9_extauc$Performance[1,], mm10_extauc$Performance[1,])
Summary.AUC$Model <- c('', 'Standard MLR', 'Ridge MLR',  'Random Forest', 'XGBoost','Neural Network', 'Support vector machine')

Summary.PI <- rbind(NA.forest, mm1_extauc$Performance[2,], mm4_extauc$Performance[2,], 
                    mm7_extauc$Performance[2,], mm8_extauc$Performance[2,],mm9_extauc$Performance[2,],mm10_extauc$Performance[2,])
Summary.PI$Model <- c('', 'Standard MLR',  'Ridge MLR', 'Random Forest', 'XGBoost','Neural Network', 'Support vector machine')

Summary.AUC
tabletext <- cbind(
  c('Model', 'Standard MLR', 'Ridge MLR',  'Random Forest', 'XGBoost','Neural Network', 'Support vector machine'),
  c('AUROC (95% CI)', 
    paste(format(round(mm1_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm1_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm1_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm4_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm4_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm4_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm7_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm7_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm7_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm8_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm8_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm8_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm9_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm9_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm9_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(mm10_extauc$Performance$AUC[1], 2), nsmall = 2), " (", format(round(mm10_extauc$Performance$LL[1], 2), nsmall = 2), " to ", format(round(mm10_extauc$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(mm1_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm1_extauc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm4_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm4_extauc$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(mm7_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm7_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm8_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm8_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm9_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm9_extauc$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(mm10_extauc$Performance$LL[2], 2), nsmall = 2), " to ", format(round(mm10_extauc$Performance$UL[2], 2), nsmall = 2), ")"))
)


tiff("Binary AUC/Summary plot.tiff", width = 31, height = 20, units = "cm", res = 300)

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
           zero=NA,
           xticks = c(0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.85, 1))
dev.off()

#####################
## Clinical utility##
#####################

#MLR
CULR1 <- DataWinBugs.imp.CM(pred=MM1_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)$`Results`

#Ridge MLR
CULR4 <- DataWinBugs.imp.CM(pred=MM4_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)$`Results`

#RF
CULR7 <- DataWinBugs.imp.CM(pred=MM7_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)$`Results`

#XGBoost
CULR8 <- DataWinBugs.imp(pred=MM8_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)$`Results`

#NN
CULR9 <- DataWinBugs.imp(pred=MM9_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)$`Results`

#SVM
CULR10 <- DataWinBugs.imp(pred=MM10_mal, outcome=outcome1, center=Center, data=testdf, imp=imp)$`Results`

# plot

tiff("NB.tiff", width = 15, height = 15, units = "cm", res = 300)
plot(datanb$threshold,datanb$`NB treatall`, type = "l", col = "gray50", lty=1, xlab = "Threshold probability of any malignancy", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.01, 0.35), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(datanb$threshold, datanb$`NB LR1`, type='l', col = "#3399CC", lty = 2, lwd = 2)
points(datanb$threshold,datanb$`NB RIDGE`, type = 'l', col = "#CCFF66", lty = 3, lwd = 3)
points(datanb$threshold, datanb$`NB RF`, type = 'l', col = "#332288", lty = 6, lwd = 3)
points(datanb$threshold, datanb$`NB XGB`, type = 'l', col = "#00FF66", lty=1,lwd = 1)
points(datanb$threshold, datanb$`NB NN`, type = 'l', col = "#FF33CC", lty=2,lwd = 2)
points(datanb$threshold, datanb$`NB SVM`, type='l', col = "#882255", lty = 3, lwd = 3)
abline(a = 0, b = 0, lty = 4, lwd = 1, col = "gray50")
legend(x=0.05, y=0.05,legend = c("Standard MLR", "Ridge MLR",  "Random Forest", "XGBoost", "Neural network", "Support vector machine","Treat all","Treat none"), ncol = 3, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c( "#3399CC", "#CCFF66",   '#332288', '#00FF66','#FF33CC',"#882255","gray50", 'gray50'), 
       lty = c(2,3,6,1,2,3,1,4), lwd = c(2,3,3,1,2,3,1,1), cex = 0.65, bty = "n")
dev.off()

##############################
### Multiclass performance ###
##############################

#Combine smallest centers into one
levels(testdf$Center)[levels(testdf$Center)=="IUK"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="MIT"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="FLI"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="NUK"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="FIT"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="CEG"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="SIT"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="UDI"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="BAI"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="LPO"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="NCI"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="PSP"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="TIT"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="MPO"] <-"Other"
levels(testdf$Center)[levels(testdf$Center)=="BCH"] <-"Other"

#########
## PDI ##
#########
#MLR
m1extpdi <- PDIimp.centre(pred=MM1pred, outcome=outcome5, center=Center,imp=imp,  data=testdf)

#Ridge
m2extpdi <- PDIimp.centre(pred=MM4pred, outcome=outcome5, center=Center,imp=imp,  data=testdf)

#RF
m5extpdi <- PDIimp.centre(pred=MM7pred, outcome=outcome5, center=Center,imp=imp,  data=testdf)

#XGBoost
m6extpdi <- PDIimp.centre(pred=pMM8, outcome=outcome5, center=Center,imp=imp,  data=testdf)

#NN
m7extpdi <- PDIimp.centre(pred=pMM9, outcome=outcome5, center=Center,imp=imp,  data=testdf)

#SVM
m8extpdi <- PDIimp.centre(pred=MM10pred, outcome=outcome5, center=Center,imp=imp,  data=testdf)


#Overall forest plot
NA.forest <- m1extpdi$Performance[1,]
NA.forest <- NA

Summary.PDI <- rbind(NA.forest, m1extpdi$Performance[1,],m2extpdi$Performance[1,], m5extpdi$Performance[1,], 
                     m6extpdi$Performance[1,],m7extpdi$Performance[1,],m8extpdi$Performance[1,])
Summary.PDI$Model <- c('',  'Standard MLR', 'Ridge MLR',   "Random forest", "XGBoost", "Neural network", "Support vector machine")

Summary.PI <- rbind(NA.forest, m1extpdi$Performance[2,], m2extpdi$Performance[2,], m5extpdi$Performance[2,], 
                    m6extpdi$Performance[2,],
                    m7extpdi$Performance[2,],m8extpdi$Performance[2,]
)
Summary.PI$Model <- c('',  'Standard MLR', 'Ridge MLR',   "Random forest", "XGBoost", "Neural network", "Support vector machine")


tabletext <- cbind(
  c('Model',  'Standard MLR', 'Ridge MLR',   "Random Forest", "XGBoost", "Neural network", "Support vector machine"),
  c('PDI (95% CI)', 
    paste(format(round(m1extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m1extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m1extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m2extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m2extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m2extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m5extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m5extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m5extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m6extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m6extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m6extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m7extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m7extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m7extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(m8extpdi$Performance$PDI[1], 2), nsmall = 2), " (", format(round(m8extpdi$Performance$LL[1], 2), nsmall = 2), " to ", format(round(m8extpdi$Performance$UL[1], 2), nsmall = 2), ")", sep = "")
  ),
  c('95% PI', 
    paste0("(", format(round(m1extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m1extpdi$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(m2extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m2extpdi$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(m5extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m5extpdi$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(m6extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m6extpdi$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(m7extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(round(m7extpdi$Performance$UL[2], 2), nsmall = 2), ")"),
    paste0("(", format(round(m8extpdi$Performance$LL[2], 2), nsmall = 2), " to ", format(0.47), ")")
  )
)


tiff("PDI/summary plot.tiff", width = 31, height = 20, units = "cm", res = 300)
forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.PDI$PDI, 3),
           lower = round(Summary.PDI$LL, 3),
           upper = round(Summary.PDI$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE),
           xlab = "PDI (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(8, "cm"),
           graph.pos = 3,
           zero=NA,
           xticks = c(0.4, 0.5, 0.6, 0.7), xlog = F, clip = c(0.4,0.7))

dev.off()

####################
### Pairwise AUC ###
####################
#Pair 1 benign vs borderline
#MLR
testdf$bvsbridge <- MM1pred[, 1] / (MM1pred[,1] + MM1pred[,2])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==2,])

#Ridge
testdf$bvsbridge <- MM4pred[, 1] / (MM4pred[,1] + MM4pred[,2])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==2,])

#RF
testdf$bvsbridge <- MM7pred[, 1] / (MM7pred[,1] + MM7pred[,2])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==2,])

#XGBoost
testdf$bvsbridge <- pMM8[, 1] / (pMM8[,1] + pMM8[,2])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==2,])

#NN
testdf$bvsbridge <- pMM9[, 1] / (pMM9[,1] + pMM9[,2])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==2,])

#SVM
testdf$bvsbridge <- MM10pred[, 1] / (MM10pred[,1] + MM10pred[,2])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==2,])

#Pair 2 benign vs stage 1
#MLR
testdf$bvsbridge <- MM1pred[, 1] / (MM1pred[,1] + MM1pred[,3])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==3,])

#Ridge
testdf$bvsbridge <- MM4pred[, 1] / (MM4pred[,1] + MM4pred[,3])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==3,])

#RF
testdf$bvsbridge <- MM7pred[, 1] / (MM7pred[,1] + MM7pred[,3])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==3,])

#XGBoost
testdf$bvsbridge <- pMM8[, 1] / (pMM8[,1] + pMM8[,3])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==3,])

#NN
testdf$bvsbridge <- pMM9[, 1] / (pMM9[,1] + pMM9[,3])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==3,])

#SVM
testdf$bvsbridge <- MM10pred[, 1] / (MM10pred[,1] + MM10pred[,3])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==3,])


#Pair 3 benign vs stage 2-4
#MLR
testdf$bvsbridge <- MM1pred[, 1] / (MM1pred[,1] + MM1pred[,4])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==4,])

#Ridge
testdf$bvsbridge <- MM4pred[, 1] / (MM4pred[,1] + MM4pred[,4])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==4,])

#RF
testdf$bvsbridge <- MM7pred[, 1] / (MM7pred[,1] + MM7pred[,4])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==4,])

#XGBoost
testdf$bvsbridge <- pMM8[, 1] / (pMM8[,1] + pMM8[,4])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==4,])

#NN
testdf$bvsbridge <- pMM9[, 1] / (pMM9[,1] + pMM9[,4])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==4,])

#SVM
testdf$bvsbridge <- MM10pred[, 1] / (MM10pred[,1] + MM10pred[,4])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==4,])


#Pair 4 benign vs secondary metastatic
#MLR
testdf$bvsbridge <- MM1pred[, 1] / (MM1pred[,1] + MM1pred[,5])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==5,])

#Ridge
testdf$bvsbridge <- MM4pred[, 1] / (MM4pred[,1] + MM4pred[,5])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==5,])

#RF
testdf$bvsbridge <- MM7pred[, 1] / (MM7pred[,1] + MM7pred[,5])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==5,])

#XGBoost
testdf$bvsbridge <- pMM8[, 1] / (pMM8[,1] + pMM8[,5])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==5,])

#NN
testdf$bvsbridge <- pMM9[, 1] / (pMM9[,1] + pMM9[,5])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==5,])

#SVM
testdf$bvsbridge <- MM10pred[, 1] / (MM10pred[,1] + MM10pred[,5])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==1|testdf$outcome5==5,])

#Pair 5
#MLR
testdf$bvsbridge <- MM1pred[, 2] / (MM1pred[,2] + MM1pred[,3])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==3,])

#Ridge
testdf$bvsbridge <- MM4pred[, 2] / (MM4pred[,2] + MM4pred[,3])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==3,])

#RF
testdf$bvsbridge <- MM7pred[, 2] / (MM7pred[,2] + MM7pred[,3])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==3,])

#XGBoost
testdf$bvsbridge <- pMM8[, 2] / (pMM8[,2] + pMM8[,3])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==3,])

#NN
testdf$bvsbridge <- pMM9[, 2] / (pMM9[,2] + pMM9[,3])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==3,])

#SVM
testdf$bvsbridge <- MM10pred[, 2] / (MM10pred[,2] + MM10pred[,3])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==3,])

#Pair 6
#MLR
testdf$bvsbridge <- MM1pred[, 2] / (MM1pred[,2] + MM1pred[,4])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==4,])

#Ridge
testdf$bvsbridge <- MM4pred[, 2] / (MM4pred[,2] + MM4pred[,4])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==4,])

#RF
testdf$bvsbridge <- MM7pred[, 2] / (MM7pred[,2] + MM7pred[,4])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==4,])

#XGBoost
testdf$bvsbridge <- pMM8[, 2] / (pMM8[,2] + pMM8[,4])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==4,])

#NN
testdf$bvsbridge <- pMM9[, 2] / (pMM9[,2] + pMM9[,4])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==4,])

#SVM
testdf$bvsbridge <- MM10pred[, 2] / (MM10pred[,2] + MM10pred[,4])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==4,])

#Pair 7
#MLR
testdf$bvsbridge <- MM1pred[, 2] / (MM1pred[,2] + MM1pred[,5])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==5,])

#Ridge
testdf$bvsbridge <- MM4pred[, 2] / (MM4pred[,2] + MM4pred[,5])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==5,])

#RF
testdf$bvsbridge <- MM7pred[, 2] / (MM7pred[,2] + MM7pred[,5])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==5,])

#XGBoost
testdf$bvsbridge <- pMM8[, 2] / (pMM8[,2] + pMM8[,5])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==5,])

#NN
testdf$bvsbridge <- pMM9[, 2] / (pMM9[,2] + pMM9[,5])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==5,])

#SVM
testdf$bvsbridge <- MM10pred[, 2] / (MM10pred[,2] + MM10pred[,5])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==2|testdf$outcome5==5,])

#Pair 8
#MLR
testdf$bvsbridge <- MM1pred[, 3] / (MM1pred[,3] + MM1pred[,4])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==4,])

#Ridge
testdf$bvsbridge <- MM4pred[, 3] / (MM4pred[,3] + MM4pred[,4])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==4,])

#RF
testdf$bvsbridge <- MM7pred[, 3] / (MM7pred[,3] + MM7pred[,4])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==4,])

#XGBoost
testdf$bvsbridge <- pMM8[, 3] / (pMM8[,3] + pMM8[,4])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==4,])

#NN
testdf$bvsbridge <- pMM9[, 3] / (pMM9[,3] + pMM9[,4])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==4,])

#SVM
testdf$bvsbridge <- MM10pred[, 3] / (MM10pred[,3] + MM10pred[,4])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==4,])

#Pair 9
#MLR
testdf$bvsbridge <- MM1pred[, 3] / (MM1pred[,3] + MM1pred[,5])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==5,])

#Ridge
testdf$bvsbridge <- MM4pred[, 3] / (MM4pred[,3] + MM4pred[,5])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==5,])

#RF
testdf$bvsbridge <- MM7pred[, 3] / (MM7pred[,3] + MM7pred[,5])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==5,])

#XGBoost
testdf$bvsbridge <- pMM8[, 3] / (pMM8[,3] + pMM8[,5])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==5,])

#NN
testdf$bvsbridge <- pMM9[, 3] / (pMM9[,3] + pMM9[,5])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==5,])

#SVM
testdf$bvsbridge <- MM10pred[, 3] / (MM10pred[,3] + MM10pred[,5])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==3|testdf$outcome5==5,])

#Pair 10
#MLR
testdf$bvsbridge <- MM1pred[, 4] / (MM1pred[,4] + MM1pred[,5])
bvsbAUC1 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==4|testdf$outcome5==5,])

#Ridge
testdf$bvsbridge <- MM4pred[, 4] / (MM4pred[,4] + MM4pred[,5])
bvsbAUC2 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==4|testdf$outcome5==5,])

#RF
testdf$bvsbridge <- MM7pred[, 4] / (MM7pred[,4] + MM7pred[,5])
bvsbAUC5 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==4|testdf$outcome5==5,])

#XGBoost
testdf$bvsbridge <- pMM8[, 4] / (pMM8[,4] + pMM8[,5])
bvsbAUC6 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==4|testdf$outcome5==5,])

#NN
testdf$bvsbridge <- pMM9[, 4] / (pMM9[,4] + pMM9[,5])
bvsbAUC7 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==4|testdf$outcome5==5,])

#SVM
testdf$bvsbridge <- MM10pred[, 4] / (MM10pred[,4] + MM10pred[,5])
bvsbAUC8 <- AUCimp.IOTA.BAYES(pred=bvsbridge, outcome=outcome5, center=Center,imp=imp, data=testdf[testdf$outcome5==4|testdf$outcome5==5,])



gg.wica2 <- dfwica %>%
  ggplot(aes(x=Outcome, y=value))+
  geom_jitter(color="grey", size=0.1, alpha=0.8)+
  #geom_point(color="grey", size=0.1, alpha=0.2)+
  geom_boxplot(outlier.shape=NA, lwd=0.6, fatten=0.6, alpha=0.2)+ #outlier.shape=NA
  #theme_ipsum()+
  xlab("")+
  ylab("Range of estimated probabilities across the models")+ theme_light()+
  ylim(0,1)+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank()) #,plot.margin = unit(c(0,0,0,0), "lines")

gg.wica2

tiff("Descriptive/boxplots of range with and without ca.tiff", width = 18, height = 20, units = "cm", res = 300)
figure <- ggarrange(gg.wica2 + rremove("ylab") , gg.woca + rremove("ylab")+rremove("xlab"), # remove axis labels from plots
                    ncol = 1, nrow = 2,
                    #align = "hv", 
                    #font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"), 
                    labels=c("A","B"), label.x=0.95, label.y=0.99)
#figure
annotate_figure(figure, left = textGrob("Range of estimated probabilities across the models", rot = 90, vjust = 1, gp = gpar(cex = 1.2)),
                bottom = textGrob("Estimated probability", gp = gpar(cex = 1.2)))
dev.off()