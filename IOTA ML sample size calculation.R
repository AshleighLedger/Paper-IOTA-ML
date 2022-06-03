#Sample size calculation IOTA ML
#########################
#Sample size calculation#
#########################
# there are 945 malignant cases out of 2489, 38% of malignancy
#step 1: sample size for O/E

#STEP (i): O/E calculation
#target confidence interval width of 1 for O/E; assuming O/E is 1 then this
#corresponds to SE(lnOE) of 0.245


#inputs, assuming an O/E =1 and width of CI is 0.1
selnoe =0.026
outcome_prop = 0.38
#calculation
sampsize_OE = (1- outcome_prop)/(outcome_prop*(selnoe^2))
sampsize_OE*outcome_prop

###########################################################################
###The second step is to consider the sample size for the calibration slope


####
# Part A, for the calibration model, specify the anticipated values of alpha and beta in the external validation population
beta0 = 0
beta1 = 1

testdf$LPMLR <- car::logit(testdf$MM1_mal)

### when specifying a separate distribution
library(fGarch)
#malignant
plot(density(testdf$LPMLR[testdf$outcome1==1]), col="red")
LPmal <- rsnorm(n = 5000000, mean = 1.5, sd = 2.4, xi = -1.1)
lines(density(LPmal), col="black")
is.na(M1$LP[M1$binaryCD==0])
#benign
plot(density(testdf$LPMLR[testdf$outcome1==0]), col="red")
LPben <- rsnorm(n = 5000000, mean = -3.3, sd = 1.2, xi = 1.1)
lines(density(LPben), col="black")


#generate dataset

df <- matrix(ncol = 2, nrow = 5000000)
colnames(df) <- c('y', 'LP')
df <- data.frame(df)
df$y[1:1900000] <- 1
df$y[1900001:5000000] <- 0
df$LP[1:1900000] <- rsnorm(n=1900000,mean = 1.5, sd = 2.4, xi = -1.1)
df$LP[1900001:5000000] <- rsnorm(n = 3100000, mean = -3.3, sd = 1.2, xi = 1.1)



LP <- df$LP



#Part D, For each participants, use their generated value of LPi from part C and the chosen values of alpha and beta from part A
#to calculate values for each ai, bi, ci

#Calculate elements of I matrix
Borenstein_00 = exp(beta0 + (beta1*LP))/((1+ exp(beta0 + (beta1*LP)))^2) #a
Borenstein_01 = LP*exp(beta0 + (beta1*LP))/((1+ exp(beta0 + (beta1*LP)))^2) #b
Borenstein_11 = LP*LP*exp(beta0 + (beta1*LP))/((1+ exp(beta0 + (beta1*LP)))^2) #c

#Part E, using all participants in the dataset, calculate the mean value of a, the mean value of b, the mean value of c
# which corresponds to I(alpha), I(alpha,beta), I(beta) respectively

Ia = mean(Borenstein_00)
Iab = mean(Borenstein_01)
Ib = mean(Borenstein_11)

#Part D, For each participants, use their generated value of LPi from part C and the chosen values of alpha and beta from part A
#to calculate values for each ai, bi, ci

#Calculate elements of I matrix
Borenstein_00 = exp(beta0 + (beta1*LP))/((1+ exp(beta0 + (beta1*LP)))^2) #a
Borenstein_01 = LP*exp(beta0 + (beta1*LP))/((1+ exp(beta0 + (beta1*LP)))^2) #b
Borenstein_11 = LP*LP*exp(beta0 + (beta1*LP))/((1+ exp(beta0 + (beta1*LP)))^2) #c

#Part E, using all participants in the dataset, calculate the mean value of a, the mean value of b, the mean value of c
# which corresponds to I(alpha), I(alpha,beta), I(beta) respectively

Ia = mean(Borenstein_00)
Iab = mean(Borenstein_01)
Ib = mean(Borenstein_11)

#for a CI width of 0.16
seslope = 0.16/(2*1.96) 
outcome_prop = 0.38
samplesize_slope = (Ia/((seslope^2)*((Ia*Ib)-(Iab^2))))

############################################################################################
###3rd step c statistics, no closed form solution so need an iterative or deductive approach
#based on the AUC of the development data, we have an AUC of 0.94. we can then expect an AUC of 0.90 for the testing data

#targeting a CI width of 0.10, the SE(C) is 0.013

outcome_prop=0.38
Cstat=0.90

#No closed form solution, need an iterative or deductive approach
#First calculate the SE for a range of N
df <- data.frame(1:10000000)
size = c(1:10000000)
seCstatq = Cstat*(1-Cstat)*(1+(((size/2)-1)*((1-Cstat)/(2-Cstat)))+
                              ((((size/2)-1)*Cstat)/(1+Cstat)))/((size^2)*outcome_prop*(1-outcome_prop))
seCstat=sqrt(seCstatq)
df$CIwidth = 2*1.96*seCstat 
df = df[df$CIwidth<=0.05,]
sampsizeCstat = min(df$X1.1e.07)

###################################
#4th step: Net Benefit calculation#



# 1% threshold
outcome_prop=0.38
sens=0.99
spec=0.128
threshold=0.01
NB=(sens*outcome_prop)-((1-spec)*(1-outcome_prop)*(threshold/(1-threshold)))

sNB=NB/outcome_prop

#Calculate NB and sNB
net_benefit=NB
standardised_net_benefit=sNB
w=((1-outcome_prop)/outcome_prop)*(threshold/(1-threshold))


#target width 0.10
sesNB=0.05/(2*1.96)
#Calculate sample size
sampsize_sNB = (1/(sesNB^2))*(((sens*(1-sens))/outcome_prop)+(w*w*spec*(1-spec)/(1-outcome_prop))+(w*w*(1-spec)*(1-spec)/(outcome_prop*(1-outcome_prop))))
sampsize_sNB



## threshold of 0.03##
outcome_prop=0.38
sens=0.95
spec=0.42
threshold=0.03
NB=(sens*outcome_prop)-((1-spec)*(1-outcome_prop)*(threshold/(1-threshold)))
sNB=NB/outcome_prop

#Calculate NB and sNB
net_benefit=NB
standardised_net_benefit=sNB
w=((1-outcome_prop)/outcome_prop)*(threshold/(1-threshold))

# target CI width for sNB of 0.1
sesNB=0.05/(2*1.96)
#Calculate sample size
sampsize_sNB = (1/(sesNB^2))*(((sens*(1-sens))/outcome_prop)+(w*w*spec*(1-spec)/(1-outcome_prop))+(w*w*(1-spec)*(1-spec)/(outcome_prop*(1-outcome_prop))))
sampsize_sNB



### threshold of 0.05

outcome_prop=0.38
sens=0.93
spec=0.61
threshold=0.05
NB=(sens*outcome_prop)-((1-spec)*(1-outcome_prop)*(threshold/(1-threshold)))
sNB=NB/outcome_prop

#Calculate NB and sNB
net_benefit=NB
standardised_net_benefit=sNB
w=((1-outcome_prop)/outcome_prop)*(threshold/(1-threshold))

# target CI width for sNB of 0.1
sesNB=0.05/(2*1.96)
#Calculate sample size
sampsize_sNB = (1/(sesNB^2))*(((sens*(1-sens))/outcome_prop)+(w*w*spec*(1-spec)/(1-outcome_prop))+(w*w*(1-spec)*(1-spec)/(outcome_prop*(1-outcome_prop))))
sampsize_sNB


### threshold of 0.10

outcome_prop=0.38
sens=0.920
spec=0.740
threshold=0.10
NB=(sens*outcome_prop)-((1-spec)*(1-outcome_prop)*(threshold/(1-threshold)))
sNB=NB/outcome_prop

#Calculate NB and sNB
net_benefit=NB
standardised_net_benefit=sNB
w=((1-outcome_prop)/outcome_prop)*(threshold/(1-threshold))

# target CI width for sNB of 0.1
sesNB=0.05/(2*1.96)
#Calculate sample size
sampsize_sNB = (1/(sesNB^2))*(((sens*(1-sens))/outcome_prop)+(w*w*spec*(1-spec)/(1-outcome_prop))+(w*w*(1-spec)*(1-spec)/(outcome_prop*(1-outcome_prop))))
sampsize_sNB

### threshold of 0.20

outcome_prop=0.38
sens=0.87
spec=0.84
threshold=0.20
NB=(sens*outcome_prop)-((1-spec)*(1-outcome_prop)*(threshold/(1-threshold)))
sNB=NB/outcome_prop

#Calculate NB and sNB
net_benefit=NB
standardised_net_benefit=sNB
w=((1-outcome_prop)/outcome_prop)*(threshold/(1-threshold))

# target CI width for sNB of 0.1
sesNB=0.05/(2*1.96)
#Calculate sample size
sampsize_sNB = (1/(sesNB^2))*(((sens*(1-sens))/outcome_prop)+(w*w*spec*(1-spec)/(1-outcome_prop))+(w*w*(1-spec)*(1-spec)/(outcome_prop*(1-outcome_prop))))
sampsize_sNB


### threshold of 0.30

outcome_prop=0.38
sens=0.80
spec=0.88
threshold=0.30
NB=(sens*outcome_prop)-((1-spec)*(1-outcome_prop)*(threshold/(1-threshold)))
sNB=NB/outcome_prop

#Calculate NB and sNB
net_benefit=NB
standardised_net_benefit=sNB
w=((1-outcome_prop)/outcome_prop)*(threshold/(1-threshold))

# target CI width for sNB of 0.1
sesNB=0.10/(2*1.96)
#Calculate sample size
sampsize_sNB = (1/(sesNB^2))*(((sens*(1-sens))/outcome_prop)+(w*w*spec*(1-spec)/(1-outcome_prop))+(w*w*(1-spec)*(1-spec)/(outcome_prop*(1-outcome_prop))))
sampsize_sNB

### threshold of 0.40

outcome_prop=0.38
sens=0.75
spec=0.91
threshold=0.40
NB=(sens*outcome_prop)-((1-spec)*(1-outcome_prop)*(threshold/(1-threshold)))
sNB=NB/outcome_prop

#Calculate NB and sNB
net_benefit=NB
standardised_net_benefit=sNB
w=((1-outcome_prop)/outcome_prop)*(threshold/(1-threshold))

# target CI width for sNB of 0.1
sesNB=0.10/(2*1.96)
#Calculate sample size
sampsize_sNB = (1/(sesNB^2))*(((sens*(1-sens))/outcome_prop)+(w*w*spec*(1-spec)/(1-outcome_prop))+(w*w*(1-spec)*(1-spec)/(outcome_prop*(1-outcome_prop))))
sampsize_sNB

### threshold of 0.50

outcome_prop=0.38
sens=0.71
spec=0.93
threshold=0.50
NB=(sens*outcome_prop)-((1-spec)*(1-outcome_prop)*(threshold/(1-threshold)))
sNB=NB/outcome_prop

#Calculate NB and sNB
net_benefit=NB
standardised_net_benefit=sNB
w=((1-outcome_prop)/outcome_prop)*(threshold/(1-threshold))

# target CI width for sNB of 0.1
sesNB=0.10/(2*1.96)
#Calculate sample size
sampsize_sNB = (1/(sesNB^2))*(((sens*(1-sens))/outcome_prop)+(w*w*spec*(1-spec)/(1-outcome_prop))+(w*w*(1-spec)*(1-spec)/(outcome_prop*(1-outcome_prop))))
sampsize_sNB
