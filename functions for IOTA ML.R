#FUNCTIONS USED 

#PDI extended from package mcca
pdi=function(y,d,k=3,...){
  num=k
  
  
  if(num==3){
    y=as.numeric(y)
    d=data.matrix(d)
    n1=which(y==1) #return the label
    n2=which(y==2)
    n3=which(y==3)
    
    #define the id
    
    pp_sum <- apply(d,1,sum)
    a <- pp_sum<0.999 | pp_sum>1.001
    b <- sum(a)
    if (b!=0){
      cat("ERROR: The input value \"d\" should be a probability matrix.")
      return(NULL)
    }
    pp=d
    
    
    pv=pp
    pv1=pv[n1,]
    pv2=pv[n2,]
    pv3=pv[n3,]
    pdi1<-0
    pdi2<-0
    pdi3<-0
    for(i in 1:length(n1)){
      pdi1=pdi1+sum(pv1[i,1]>pv2[,1])*sum(pv1[i,1]>pv3[,1])
    }
    for(i in 1:length(n2)){
      pdi2=pdi2+sum(pv2[i,2]>pv1[,2])*sum(pv2[i,2]>pv3[,2])
    }
    for(i in 1:length(n3)){
      pdi3=pdi3+sum(pv3[i,3]>pv1[,3])*sum(pv3[i,3]>pv2[,3])
    }
    pdi<-(pdi1+pdi2+pdi3)/(3*length(n1)*length(n2)*length(n3))
    return(pdi)
    
  }else if(num==4){
    y=as.numeric(y)
    d=data.matrix(d)
    
    n1=which(y==1) #return the label
    n2=which(y==2)
    n3=which(y==3)
    n4=which(y==4)
    
    #define the id
    
    pp_sum <- apply(d,1,sum)
    a <- pp_sum<0.999 | pp_sum>1.001
    b <- sum(a)
    if (b!=0){
      cat("ERROR: The input value \"d\" should be a probability matrix.")
      return(NULL)
    }
    pp=d
    
    
    pv=pp
    
    pv1=pv[n1,]
    pv2=pv[n2,]
    pv3=pv[n3,]
    pv4=pv[n4,]
    pdi1<-0
    pdi2<-0
    pdi3<-0
    pdi4<-0
    for(i in 1:length(n1)){
      pdi1=pdi1+sum(pv1[i,1]>pv2[,1])*sum(pv1[i,1]>pv3[,1])*sum(pv1[i,1]>pv4[,1])
    }
    for(i in 1:length(n2)){
      pdi2=pdi2+sum(pv2[i,2]>pv1[,2])*sum(pv2[i,2]>pv3[,2])*sum(pv2[i,2]>pv4[,2])
    }
    for(i in 1:length(n3)){
      pdi3=pdi3+sum(pv3[i,3]>pv1[,3])*sum(pv3[i,3]>pv2[,3])*sum(pv3[i,3]>pv4[,3])
    }
    for(i in 1:length(n4)){
      pdi4=pdi4+sum(pv4[i,4]>pv1[,4])*sum(pv4[i,4]>pv2[,4])*sum(pv4[i,4]>pv3[,4])
    }
    pdi<-(pdi1+pdi2+pdi3+pdi4)/(4*length(n1)*length(n2)*length(n3)*length(n4))
    return(pdi)
    
  }else if(num==2){
    
    y=as.numeric(y)
    d=data.matrix(d)
    n1=which(y==1) #return the label
    n2=which(y==2)
    
    
    #define the id
    
    pp_sum <- apply(d,1,sum)
    a <- pp_sum<0.999 | pp_sum>1.001
    b <- sum(a)
    if (b!=0){
      cat("ERROR: The input value \"d\" should be a probability matrix.")
      return(NULL)
    }
    pp=d
    
    
    pv=pp
    pv1=pv[n1,]
    pv2=pv[n2,]
    
    pdi1<-0
    pdi2<-0
    
    for(i in 1:length(n1)){
      pdi1=pdi1+sum(pv1[i,1]>pv2[,1])+sum(pv1[i,1]==pv2[,1])
    }
    for(i in 1:length(n2)){
      pdi2=pdi2+sum(pv2[i,2]>pv1[,2])
    }
    
    pdi<-(pdi1+pdi2)/(2*length(n1)*length(n2))
    return(pdi)
    
  }else if(num==5){
    y <- as.numeric(y)
    d <- data.matrix(d)
    n1 <- which(y==1)
    n2 <- which(y==2)
    n3 <- as.integer(which(y==3))
    n4 <- which(y==4)
    n5 <- which(y==5)
    pp_sum <- apply(d,1,sum)
    a <- pp_sum<0.999 | pp_sum>1.001
    b <- sum(a)
    if (b!=0){
      cat("ERROR: The input value \"d\" should be a probability matrix.")
      return(NULL)
    }
    pp=d
    pv=pp
    
    pv1=pv[n1,]
    pv2=pv[n2,]
    pv3=pv[n3,]
    pv4=pv[n4,]
    pv5=pv[n5,]
    pdi1<-0
    pdi2<-0
    pdi3<-0
    pdi4<-0
    pdi5<-0
    for(i in 1:length(n1)){
      pdi1=pdi1+as.double(sum(pv1[i,1]>pv2[,1]))*as.double(sum(pv1[i,1]>pv3[,1]))*as.double(sum(pv1[i,1]>pv4[,1]))*as.double(sum(pv1[i,1]>pv5[,1]))
    }
    for(i in 1:length(n2)){
      pdi2=pdi2+as.double(sum(pv2[i,2]>pv1[,2]))*as.double(sum(pv2[i,2]>pv3[,2]))*as.double(sum(pv2[i,2]>pv4[,2]))*as.double(sum(pv2[i,2]>pv5[,2]))
    }
    for(i in 1:length(n3)){
      pdi3=pdi3+as.double(sum(pv3[i,3]>pv1[,3]))*as.double(sum(pv3[i,3]>pv2[,3]))*as.double(sum(pv3[i,3]>pv4[,3]))*as.double(sum(pv3[i,3]>pv5[,3]))
    }
    for(i in 1:length(n4)){
      pdi4=pdi4+as.double(sum(pv4[i,4]>pv1[,4]))*as.double(sum(pv4[i,4]>pv2[,4]))*as.double(sum(pv4[i,4]>pv3[,4]))*as.double(sum(pv4[i,4]>pv5[,4]))
    }
    for(i in 1:length(n5)){
      pdi5=pdi5+as.double(sum(pv5[i,5]>pv1[,5]))*as.double(sum(pv5[i,5]>pv2[,5]))*as.double(sum(pv5[i,5]>pv3[,5]))*as.double(sum(pv5[i,5]>pv4[,5]))
    }
    pdi<-(pdi1+pdi2+pdi3+pdi4+pdi5)/(5*length(n1)*length(n2)*length(n3)*length(n4)*length(n5))
    return(pdi)
  }
  
}


ests <- function (y, d, acc="hum",level=0.95,method = "multinom", k = 3,B=50,balance=FALSE, ...) {
  
  series=numeric()
  
  if (acc=="hum"){
    
    if (balance==FALSE){
      for (b in 1:B){
        nn <- length(y)
        id <- sample(1:nn,nn,replace = T)
        #id <- unique(id)
        while (length(unique(y[id]))<k){
          id <- sample(1:nn,nn,replace = T)
          #id <- unique(id)
        }
        while (min(table(y[id]))<2){
          id <- sample(1:nn,nn,replace = T)
          #id <- unique(id)
        }
        if (class(d)=="numeric"){
          series[b] <- hum(y=y[id],d=d[id],method=method,k=k,...)
        }else {
          series[b] <- hum(y=y[id],d=d[id,],method=method,k=k,...)
        }
      }
    }
    if (balance==TRUE){
      for (b in 1:B){
        id <- unlist(caret::createResample(y, times = 1))
        if (class(d)=="numeric"){
          series[b] <- hum(y=y[id],d=d[id],method=method,k=k,...)
        }else {
          series[b] <- hum(y=y[id],d=d[id,],method=method,k=k,...)
        }
      }
    }
    
    series.sort <- sort(series)
    return(list(value=hum(y=y,d=d,method=method,k=k,...),
                se=sd(series),
                interval=c(series.sort[ifelse(B*(0.5-level/2)<1,1,B*(0.5-level/2))],series.sort[B*(0.5+level/2)])))
  }
  if (acc=="pdi"){
    if (balance==FALSE){
      for (b in 1:B){
        nn <- length(y)
        id <- sample(1:nn,nn,replace = T)
        #id <- unique(id)
        while (length(unique(y[id]))<k){
          id <- sample(1:nn,nn,replace = T)
          #id <- unique(id)
        }
        while (min(table(y[id]))<2){
          id <- sample(1:nn,nn,replace = T)
          #id <- unique(id)
        }
        if (class(d)=="numeric"){
          series[b] <- pdi(y=y[id],d=d[id],method=method,k=k,...)
        }else {
          series[b] <- pdi(y=y[id],d=d[id,],method=method,k=k,...)
        }
      }
    }
    if (balance==TRUE){
      for (b in 1:B){
        id <- unlist(caret::createResample(y, times = 1))
        if (class(d)=="numeric"){
          series[b] <- pdi(y=y[id],d=d[id],method=method,k=k,...)
        }else {
          series[b] <- pdi(y=y[id],d=d[id,],method=method,k=k,...)
        }
      }
    }
    series.sort <- sort(series)
    return(list(value=pdi(y=y,d=d,method=method,k=k,...),
                se=sd(series),
                interval=c(series.sort[ifelse(B*(0.5-level/2)<1,1,B*(0.5-level/2))],series.sort[B*(0.5+level/2)])))
  }
  if (acc=="ccp"){
    if (balance==FALSE){
      for (b in 1:B){
        nn <- length(y)
        id <- sample(1:nn,nn,replace = T)
        #id <- unique(id)
        while (length(unique(y[id]))<k){
          id <- sample(1:nn,nn,replace = T)
          #id <- unique(id)
        }
        while (min(table(y[id]))<2){
          id <- sample(1:nn,nn,replace = T)
          #id <- unique(id)
        }
        if (class(d)=="numeric"){
          series[b] <- ccp(y=y[id],d=d[id],method=method,k=k,...)
        }else {
          series[b] <- ccp(y=y[id],d=d[id,],method=method,k=k,...)
        }
      }
    }
    if (balance==TRUE){
      for (b in 1:B){
        id <- unlist(caret::createResample(y, times = 1))
        if (class(d)=="numeric"){
          series[b] <- ccp(y=y[id],d=d[id],method=method,k=k,...)
        }else {
          series[b] <- ccp(y=y[id],d=d[id,],method=method,k=k,...)
        }
      }
    }
    series.sort <- sort(series)
    return(list(value=ccp(y=y,d=d,method=method,k=k,...),
                se=sd(series),
                interval=c(series.sort[ifelse(B*(0.5-level/2)<1,1,B*(0.5-level/2))],series.sort[B*(0.5+level/2)])))
  }
  
  if (acc=="rsq"){
    if (balance==FALSE){
      for (b in 1:B){
        nn <- length(y)
        id <- sample(1:nn,nn,replace = T)
        #id <- unique(id)
        while (length(unique(y[id]))<k){
          id <- sample(1:nn,nn,replace = T)
          #id <- unique(id)
        }
        while (min(table(y[id]))<2){
          id <- sample(1:nn,nn,replace = T)
          #id <- unique(id)
        }
        if (class(d)=="numeric"){
          series[b] <- rsq(y=y[id],d=d[id],method=method,k=k,...)
        }else {
          series[b] <- rsq(y=y[id],d=d[id,],method=method,k=k,...)
        }
      }
    }
    if (balance==TRUE){
      for (b in 1:B){
        id <- unlist(caret::createResample(y, times = 1))
        if (class(d)=="numeric"){
          series[b] <- rsq(y=y[id],d=d[id],method=method,k=k,...)
        }else {
          series[b] <- rsq(y=y[id],d=d[id,],method=method,k=k,...)
        }
      }
    }
    series.sort <- sort(series)
    return(list(value=rsq(y=y,d=d,method=method,k=k,...),
                se=sd(series),
                interval=c(series.sort[ifelse(B*(0.5-level/2)<1,1,B*(0.5-level/2))],series.sort[B*(0.5+level/2)])))
  }
  
}



######
AUCimp.IOTA <- function(pred, outcome, center, imp, data, method.MA = "REML", titleGraph = "AUC per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  
  AUCimp <- list()
  for(i in 1:NRimp){
    AUCimp[[i]] <- list()
  }
  
  PrevalenceOverall <- matrix(nrow = NRimp, ncol = 1)
  
  # AUC per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    AUCcenter <- matrix(ncol = 6, nrow = length(centers))
    colnames(AUCcenter) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
    AUCcenter <- data.frame(AUCcenter)
    AUCcenter$Center <- centers
    
    for(i in seq_along(centers)){
      AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1 & Df$imp == j], Df$p[Df$center == centers[i] & Df$y == 0 & Df$imp == j], method = "pepe")
      AUCcenter[i, 2]   <- nrow(Df[Df$center == centers[i] & Df$imp == j,])
      AUCcenter[i, 3]   <- round(nrow(Df[Df$y == 1 & Df$center == centers[i] & Df$imp == j,])/nrow(Df[Df$center == centers[i] & Df$imp == j,])*100)
      
      ## Additional part for AUCs of 1
      if(AUCcenter[i, 4] == 1){
        AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1 & Df$imp == j], Df$p[Df$center == centers[i] & Df$y == 0 & Df$imp == j], method = "newcombe") # Newcombe ipv pepe
      } else{
        AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1 & Df$imp == j], Df$p[Df$center == centers[i] & Df$y == 0 & Df$imp == j], method = "pepe")
      }
      
      if(AUCcenter$AUC[i] != 1){
        AUCcenter$logit.AUC[i] <- logit(AUCcenter$AUC[i])
        AUCcenter$logit.se[i]  <- (logit(AUCcenter$AUC[i]) - logit(AUCcenter$LL[i]))/1.96
        AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
      } else{
        AUCcenter$logit.AUC[i] <- logit(0.999)
        AUCcenter$logit.se[i]  <- (logit(0.999) - logit(AUCcenter$LL[i]))/1.96
        AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
      }
      
    }
    AUCcenter <- AUCcenter[order(AUCcenter$SampleSize, decreasing = TRUE),]
    
    AUCimp[[j]] <- AUCcenter
    
    PrevalenceOverall[j] <- round(nrow(Df[Df$y == 1 & Df$imp == j,])/nrow(Df[Df$imp == j,])*100)
  }
  
  AUCimpLong <- rbindlist(AUCimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  AUCcombined <- matrix(ncol = 6, nrow = length(centers))
  colnames(AUCcombined) <- c('Center', 'SampleSize', 'Prevalence', 'logit.AUC', 'logit.LL', 'logit.UL')
  AUCcombined <- data.frame(AUCcombined)
  AUCcombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    AUCcombined$SampleSize[i] <- unique(AUCimpLong$SampleSize[AUCimpLong$Center == centers[i]])
    AUCcombined$Prevalence[i] <- round(mean(AUCimpLong$Prevalence[AUCimpLong$Center == centers[i]]))
    AUCcombined[i, 4] <- mean(AUCimpLong$logit.AUC[AUCimpLong$Center == centers[i]])
    WithinVar <- mean(AUCimpLong$logit.var[AUCimpLong$Center == centers[i]])
    BetweenVar <- var(AUCimpLong$logit.AUC[AUCimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    AUCcombined$PooledSE[i] <- sqrt(PooledVar)
    AUCcombined$logit.LL[i] <- AUCcombined$logit.AUC[i] - 1.96*AUCcombined$PooledSE[i]
    AUCcombined$logit.UL[i] <- AUCcombined$logit.AUC[i] + 1.96*AUCcombined$PooledSE[i]
  }
  
  AUCcombined$AUC <- inv.logit(AUCcombined$logit.AUC)
  AUCcombined$LL <- inv.logit(AUCcombined$logit.LL)
  AUCcombined$UL <- inv.logit(AUCcombined$logit.UL)
  AUCcombined <- AUCcombined[order(AUCcombined$SampleSize, decreasing = TRUE),]
  
  AUCoverall <- matrix(nrow = 2, ncol = 6)
  colnames(AUCoverall) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCoverall <- data.frame(AUCoverall)
  AUCoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  AUCoverall$SampleSize <- nrow(Df[Df$imp == 1,])
  AUCoverall$Prevalence <- round(mean(PrevalenceOverall))
  
  # Meta-analyse voor overall estimate
  fit.RE = rma.uni(AUCcombined$logit.AUC, sei = AUCcombined$PooledSE, method = method.MA)
  PI = predict(fit.RE, transf = transf.ilogit)
  
  AUCoverall$AUC[1] <- inv.logit(coef(fit.RE))
  AUCoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  AUCoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  AUCoverall$AUC[2] <- PI$pred
  AUCoverall$LL[2] <- PI$cr.lb
  AUCoverall$UL[2] <- PI$cr.ub
  AUCoverall
  
  NAforest <- matrix(nrow = 1, ncol = 6)
  colnames(NAforest) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  AUC <- rbind(AUCcombined[, c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')], NAforest, NAforest, AUCoverall)
  
  # Layout for forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(AUC)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- AUC$Center[i]
  }
  RRauc <- 1:nrobs
  for(i in 1:nrobs){
    RRauc[i] <- paste(format(round(AUC$AUC[i], 2), nsmall = 2), " (", format(round(AUC$LL[i], 2), nsmall = 2), " to ", format(round(AUC$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- AUC$SampleSize[i]
  }
  
  Labels <- c('Centre', 'AUROC (95% CI)', 'N') #(prev)')
  Combined <- cbind(RRcenter, RRauc, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  AUCna <- rbind(NAforest, NAforest, AUC)
  
  return(structure(list(Performance = AUCoverall, ModelFit = fit.RE, AUCcenters = AUCcombined[, c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')], IncludedCenters = centers, dataPlot = AUCna, Plot = Tabletext, data = Df)))
  
}

rubinrules <- function (estimate, variance, data, alpha = 0.05, exp = F, expit = F) 
{
  if (missing(data)) {
    Est = estimate
    Var = variance
  }
  else {
    if (!is.data.frame(data)) 
      stop("Data must be of type dataframe")
    Argz = as.list(match.call())[-1]
    Est = eval(Argz$estimate, data)
    Var = eval(Argz$variance, data)
  }
  AvgEst = mean(Est)
  AvgVar = mean(Var)
  Nimp = length(Est)
  Bvari = crossprod(Est - AvgEst)
  Bvar = Bvari/(Nimp - 1)
  TotalVar = AvgVar + Bvar * (1 + 1/Nimp)
  TotalSE = sqrt(TotalVar)
  Tdf = (Nimp - 1) * (1 + (Bvar * (1 + 1/Nimp))/AvgVar)^2
  Tcv = qt(1 - alpha/2, Tdf)
  LCL = AvgEst - Tcv * TotalSE
  UCL = AvgEst + Tcv * TotalSE
  if (exp) {
    ExpAvgEst = exp(AvgEst)
    ExpLCL = exp(LCL)
    ExpUCL = exp(UCL)
  }
  if (expit) {
    ExpitAvgEst = inv.logit(AvgEst)
    ExpitLCL = inv.logit(LCL)
    ExpitUCL = inv.logit(UCL)
  }
  Results = list(`Mean estimate` = AvgEst, `Lower confidence limit` = LCL, 
                 `Upper confidence limit` = UCL, `Within variance` = AvgVar, 
                 `Between variance` = Bvar, `Total variance` = TotalVar, 
                 `Total SE` = TotalSE, `Degrees of freedom, t-distribution` = Tdf, 
                 `Quantile of Student's t-distribution` = Tcv)
  if (exp) 
    Results$Exponent = c(`Exponent mean estimate` = ExpAvgEst, 
                         `Exponent LCL` = ExpLCL, `Exponent UCL` = ExpUCL)
  if (expit) 
    Results$Expit = c(`Expit mean estimate` = ExpitAvgEst, 
                      `Expit LCL` = ExpitLCL, `Expit UCL` = ExpitUCL)
  return(Results)
}

ListToVector <-function (x){unname(unlist(x))}


## Probabilities (LR2, SRrisks, ADNEX): Multiple imputed data
RE.ValProbImp <- function (p, y, center, imputation.id, patientid, data, LogCal = T, 
                           flexible = F, CL = c("none", "CI", "PI"), CalibrLines = c("overall", 
                                                                                     "centers", "both"), dostats = T, statloc = c(0, 0.85), 
                           legendloc = c(0.5, 0.27), roundstats = 2, cex = 0.75, cex.leg = 0.75, 
                           ncol.leg = 1, lty.overall = 1, lwd.overall = 2, col.overall = "red", 
                           RMprompt = F, RmNonConv = F, title = "Calibration curve", xlab = "Estimated risk of malignancy", 
                           ylab = "Observed proportion of malignancy", xlim = c(0, 1), ylim = c(0, 
                                                                                                1), d0lab = "0", d1lab = "1", cex.d01 = 0.7, dist.label = 0.04, 
                           line.bins = -0.05, dist.label2 = 0.03, las = 1, length.seg = 1, 
                           y.intersp = 1, lty.ideal = 1, col.ideal = "black", lwd.ideal = 1.75, 
                           lty.centers = NULL, lwd.centers = NULL, col.centers = NULL, 
                           Parallel = F, NrCores = detectCores() - 1, alpha = 0.05, 
                           EpsGrad = 0.001, fNonConv = c("warning", "stop"), Controlglmer = glmerControl(optimizer = "bobyqa"), 
                           AUCmeth = c("RE.auc.imp", "RE.auc.imp2"), ...) 
{
  Argz = as.list(match.call())[-1]
  p = eval(Argz$p, data)
  y = eval(Argz$y, data)
  center = eval(Argz$center, data)
  id = eval(Argz$patientid, data)
  LP = logit(p)
  CL = match.arg(CL)
  CalibrLines = match.arg(CalibrLines)
  AUCmeth = match.arg(AUCmeth)
  if (is.factor(center)) 
    center = as.character(center)
  if (length(unique(center)) == 1) 
    stop("There is only one center, hence a random effects model should not be used here.")
  if (length(unique(center)) < 5) 
    warning("There are less than 5 centers, consider using a different method.", 
            immediate. = T)
  if (LogCal & flexible) 
    stop("LogCal and flexible cannot both be true.")
  if (missing(imputation.id)) 
    stop("Specify the variable indicating the number of the imputation!")
  imp = eval(Argz$imputation.id, data)
  if (!is.numeric(imp)) 
    imp = as.numeric(imp)
  NrImp = length(unique(imp))
  Df = data.frame(y = y, p = p, LP = logit(p), center = center, 
                  imp = imp, id = id)
  Df = Df[with(Df, order(center, p)), ]
  RmCenter = dlply(Df[Df$imp == 1, ], .(center), function(x) if (sum(x$y == 
                                                                     1) < 10 | sum(x$y == 0) < 10) 
    unique(x$center)
    else NULL)
  RmCenter = ListToVector(RmCenter)
  if (length(RmCenter) != 0) {
    if (RMprompt) {
      Cmmnd = readline(paste0("The center(s) ", paste0(RmCenter, 
                                                       collapse = ", "), " have less than 10 (non-)events and these will be", 
                              " removed. Do you want to continue? (y/n)   "))
      if (Cmmnd == "n") 
        stop("Function was stopped by the user.")
    }
    Df = Df[!(Df$center %in% RmCenter), ]
  }
  IncludedCenters = unique(Df$center)
  centers = IncludedCenters
  nCenters = length(IncludedCenters)
  if (!Parallel & length(unique(Df$imp)) > 50) {
    cat("\n\nParallel computing can be used to get the results faster.\n\n")
  }
  else if (Parallel & length(unique(Df$imp)) < 10) {
    cat(paste("\n\nParallel computing is only recommended when there is a large number of imputed datasets.", 
              " When used with a small number of imputed datasets, this may be slower.\n\n"))
  }
  if (CalibrLines != "overall") {
    if (all(sapply(list(lty.centers, lwd.centers, col.centers), 
                   is.null))) {
      lty.centers = rep(1, nCenters)
      lwd.centers = rep(1, nCenters)
      col.centers = seq_len(nCenters)
    }
    else if (sum(!sapply(list(lty.centers, lwd.centers, col.centers), 
                         is.null)) != 3) {
      stop("If you specify one of the arguments, please also specify the others.")
    }
    else {
      if (any(sapply(list(lty.centers, lwd.centers, col.centers), 
                     function(x) length(x) < nCenters))) 
        stop("The vector length of lty.centers, lwd.centers and col.centers is less than the number of centers.")
      FixL = function(x) x[1:nCenters]
      lty.centers = FixL(lty.centers)
      col.centers = FixL(col.centers)
      lwd.centers = FixL(lwd.centers)
    }
  }
  RubinRes <- function(x) {
    do.call("rbind", lapply(x, function(x) {
      tmp = rubinrules(Est, Var, x, alpha = alpha)[c("Mean estimate", 
                                                            "Total variance")]
      tmp = do.call("cbind.data.frame", tmp)
      names(tmp) = c("Est", "Var")
      tmp
    }))
  }
  RubinCL <- function(x) {
    do.call("rbind", lapply(x, function(x) {
      tmp = rubinrules(Est, Var, x, alpha = alpha, expit = T)$Expit
      names(tmp) = c("Est", "LCL", "UCL")
      tmp
    }))
  }
  NewX = data.frame(LP = seq(min(Df$LP), max(Df$LP), length = 500))
  PerC = data.frame(LP = rep(NewX$LP, length(IncludedCenters)), 
                    center = sort(rep(IncludedCenters, 500)))
  ArgzCalc = list(X = 1:NrImp, FUN = function(i) {
    tmp = Df[Df$imp == i, ]
    LogMM = glmer(y ~ LP + (LP | center), data = tmp, family = "binomial", 
                  control = Controlglmer)
    LogMM2 = glmer(y ~ 1 + (1 | center), data = tmp, family = "binomial", 
                   offset = LP, control = Controlglmer)
    X = if (CL != "PI") {
      as.matrix(cbind(1, NewX$LP))
    } else {
      model.matrix(LogMM)[order(tmp$p), ]
    }
    FE = fixef(LogMM)
    VarFE = diag(as.matrix(vcov.merMod(LogMM)))
    Est = X %*% FE
    Var = diag(X %*% tcrossprod(as.matrix(vcov.merMod(LogMM)), 
                                X))
    FE_int = fixef(LogMM2)
    VarFE_int = diag(as.matrix(vcov.merMod(LogMM2)))
    Results = list(B = list(Interc = cbind.data.frame(Est = FE[1], 
                                                      Var = VarFE[1]), Slope = cbind.data.frame(Est = FE[2], 
                                                                                                Var = VarFE[2]), Interc2 = cbind.data.frame(Est = FE_int[1], 
                                                                                                                                            Var = VarFE_int[1])), LP = cbind.data.frame(Est = Est, 
                                                                                                                                                                                        Var = Var), ModelFit = LogMM)
    if (CalibrLines != "overall") {
      Results$EstC = predict(LogMM, newdata = PerC, re.form = ~(LP | 
                                                                  center), allow.new.levels = T, type = "response")
    }
    return(Results)
  }, simplify = F, USE.NAMES = T)
  cat("\n\nComputing overall/center-specific calibration curve(s)...\n\n")
  if (Parallel) {
    cl <- makeCluster(NrCores)
    clusterEvalQ(cl, library("REMA"))
    clusterExport(cl, c("Df", "NrImp", "NewX", "CL", "CalibrLines", 
                        "Controlglmer", "EpsGrad", "fNonConv"), envir = environment())
    MethCalc = "parSapply"
    ArgzCalc$cl = cl
    on.exit(stopCluster(cl))
  }
  else {
    MethCalc = "sapply"
  }
  ResultsCal = do.call(MethCalc, args = ArgzCalc)
  MF = lapply(ResultsCal, "[[", "ModelFit")
  Conv = unlist(lapply(MF, function(x) length(x@optinfo$conv$lme4) == 
                         0))
  if (RmNonConv) 
    ResultsCal = ResultsCal[Conv]
  ResultsB = lapply(ResultsCal, "[[", "B")
  ResultsBint = do.call("rbind", lapply(ResultsB, "[[", "Interc"))
  ResultsBint2 = do.call("rbind", lapply(ResultsB, "[[", "Interc2"))
  ResultsBslo = do.call("rbind", lapply(ResultsB, "[[", "Slope"))
  ResultsBRR = as.vector(RubinRes(list(ResultsBint, ResultsBslo))[, 
                                                                  1])
  SumSlo = RubinRes(list(ResultsBslo))
  SumSlo$LL <- SumSlo$Est - 1.96 * sqrt(SumSlo$Var)
  SumSlo$UL <- SumSlo$Est + 1.96 * sqrt(SumSlo$Var)
  
  SumInt = RubinRes(list(ResultsBint2))
  SumInt$LL <- SumInt$Est - 1.96 * sqrt(SumInt$Var)
  SumInt$UL <- SumInt$Est + 1.96 * sqrt(SumInt$Var)
  
  ResultsLP = lapply(ResultsCal, "[[", "LP")
  Niter = if (CL != "PI") {
    nrow(NewX)
  }
  else {
    nrow(Df[Df$imp == 1, ])
  }
  AllPred = sapply(1:Niter, function(i) {
    do.call("rbind", lapply(ResultsLP, function(x) x[i, ]))
  }, simplify = F)
  ResultsCalRR = RubinCL(AllPred)
  if (CalibrLines != "overall") {
    ResultsC = lapply(ResultsCal, "[[", "EstC")
    ResultsC = rowMeans(do.call("cbind", ResultsC))
    ResultsC = cbind.data.frame(x = inv.logit(PerC$LP), y = ResultsC, 
                                center = PerC$center)
    ResultsC = split(ResultsC, ResultsC$center)
  }
  
  if(CalibrLines != "overall"){
    par(mar = c(5,5,1,13), xpd=TRUE, pty = 's') 
  } else{
    par(pty = 's')
  }
  
  plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
       ylab = ylab, las = las, main = title, cex.lab = 1, cex.axis = 1, ...) 
  clip(0, 1, 0, 1)
  abline(0, 1, lty = lty.ideal, col = col.ideal, lwd = lwd.ideal)
  do.call("clip", as.list(par()$usr))
  lt = lty.ideal
  lw.d = lwd.ideal
  all.col = col.ideal
  leg = "Ideal"
  if (CalibrLines == "centers") {
    for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                              i)], lty = lty.centers[which(IncludedCenters == i)], 
                                     lwd = lwd.centers[which(IncludedCenters == i)])
    lt = c(lt, lty.centers)
    lw.d = c(lw.d, lwd.centers)
    all.col = c(all.col, col.centers)
    leg = c(leg, as.character(IncludedCenters))
  }
  else {
    p = if (CL == "PI") {
      ddply(Df, .(id), function(x) mean(x$p))$V1
    }
    else {
      inv.logit(NewX$LP)
    }
    X = cbind(1, logit(p))
    FE = ResultsBRR
    y = Df$y[Df$imp == 1]
    if (!flexible) {
      OverallCal = inv.logit(X[order(p), ] %*% FE)
      p = p[order(p)]
      lines(p, OverallCal, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
      if (CL != "none") {
        lines(p, ResultsCalRR[, 2], lty = 2, col = col.overall)
        lines(p, ResultsCalRR[, 3], lty = 2, col = col.overall)
      }
    }
    else {
      Lfit = loess(y ~ p, Df)
      x = Lfit$x
      y = Lfit$fitted
      y = y[order(x)]
      x = x[order(x)]
      lines(x, y, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
    }
    lt = c(lt, lty.overall)
    lw.d = c(lw.d, lwd.overall)
    all.col = c(all.col, col.overall)
    leg = c(leg, "Overall")
    if (CalibrLines == "both") {
      for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                                i)], lty = lty.centers[which(IncludedCenters == 
                                                                                                               i)], lwd = lwd.centers[which(IncludedCenters == 
                                                                                                                                              i)])
      lt = c(lt, lty.centers)
      lw.d = c(lw.d, lwd.centers)
      all.col = c(all.col, col.centers)
      leg = c(leg, as.character(IncludedCenters))
    }
  }
  lp = legendloc
  lp = list(x = lp[1], y = lp[2])
  legend("topright", leg, lty = lt, cex = 1, bty = "n", lwd = lw.d, 
         col = all.col, y.intersp = y.intersp, ncol = ncol.leg, inset = c(-.58, 0), xpd=NA) 
  if (dostats) {
    if (CalibrLines == "centers"){
      stats.2 <- paste("")
      text(statloc[1], statloc[2], stats.2, pos = 4, cex = cex)
    } else{
      stats.2 <- matrix(ncol = 2, nrow = 2)
      colnames(stats.2) <- c("", "Estimate (95% CI)")
      stats.2[1, ] <- c("Intercept", paste0(format(round(SumInt$Est, 2), nsmall = 2), " (", format(round(SumInt$LL, 2), nsmall = 2), " to ", format(round(SumInt$UL, 2), nsmall = 2), ")"))
      stats.2[2, ] <- c("Slope", paste0(format(round(SumSlo$Est, 2), nsmall = 2), " (", format(round(SumSlo$LL, 2), nsmall = 2), " to ", format(round(SumSlo$UL, 2), nsmall = 2), ")"))
      
      addtable2plot(x = statloc[1], y = statloc[2], table = stats.2, display.colnames = TRUE, cex = 0.75)
      
    }
  }
  
  if (CalibrLines != "centers") {
    x = ddply(Df, .(id), function(x) mean(x$p))$V1
    y = Df$y[Df$imp == 1]
    bins <- seq(0, min(1, max(xlim)), length = 101)
    x <- x[x >= 0 & x <= 1]
    f0 = table(cut(x[y == 0], bins))
    f1 = table(cut(x[y == 1], bins))
    j0 = f0 > 0
    j1 = f1 > 0
    bins0 = (bins[-101])[j0]
    bins1 = (bins[-101])[j1]
    f0 = f0[j0]
    f1 = f1[j1]
    maxf = max(f0, f1)
    f0 = (0.1 * f0)/maxf
    f1 = (0.1 * f1)/maxf
    segments(bins1, line.bins, bins1, length.seg * f1 + line.bins)
    segments(bins0, line.bins, bins0, length.seg * -f0 + 
               line.bins)
    lines(c(min(bins0, bins1) - 0.01, max(bins0, bins1) + 
              0.01), c(line.bins, line.bins))
    text(max(bins0, bins1) + dist.label, line.bins + dist.label2, 
         d1lab, cex = cex.d01)
    text(max(bins0, bins1) + dist.label, line.bins - dist.label2, 
         d0lab, cex = cex.d01)
  }
  
  AllResults = structure(list( 
    included = unique(Df$center), ConfLevel = 1 - alpha, 
    ResultsBRR = ResultsBRR, ResultsCalRR = ResultsCalRR, 
    PlotArgz = list(Mtext = list(stats = stats.2), Plot = list(x = lp, 
                                                               leg = leg, cex.leg = cex.leg, lwd = lw.d, col = all.col, 
                                                               y.intersp = y.intersp, ncol = ncol.leg)), call = Argz, OneInt = ResultsBint2, OneSlo = ResultsBslo, 
    One = ResultsB, SumSlo = SumSlo, SumInt = SumInt),
    class = "RE_ValProbImp")
  if (CalibrLines == "centers") {
    AllResults$PerCenter = ResultsC
  }
  else {
    AllResults$Plot = cbind.data.frame(x = p, y = OverallCal)
    if (CalibrLines == "both") 
      AllResults$PerCenter = ResultsC
  }
  return(AllResults)
}


AUC.imp <- function(pred, outcome, imp, data){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, imp = imp, stringsAsFactors = F)
  
  NRimp <- length(unique(Df$imp))
  
  AUCimp <- matrix(ncol = 3, nrow = NRimp)
  colnames(AUCimp) <- c('AUC', 'LL', 'UL')
  AUCimp <- data.frame(AUCimp)
  
  # AUC per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    AUC <- auc.nonpara.mw(Df$p[Df$y == 5 & Df$imp == j], Df$p[Df$y == 4 & Df$imp == j], method = "pepe")
    
    if(AUC[1] < 0.50){
      AUC2 <- auc.nonpara.mw(Df$p[Df$y == 4 & Df$imp == j], Df$p[Df$y == 5 & Df$imp == j], method = "pepe")
      AUCimp[j, 1:3] <- AUC2 # Zet output in tabel
    } 
    
    ## Additional part for AUCs of 1
    else{
      if(AUC[1] == 1){
        AUC3 <- auc.nonpara.mw(Df$p[Df$y == 5 & Df$imp == j], Df$p[Df$y == 4 & Df$imp == j], method = "newcombe") # Newcombe ipv pepe
        AUCimp[j, 1:3] <- AUC3 # Zet output in tabel
      } else{
        AUCimp[j, 1:3] <- AUC # Zet output in tabel
      }
    }
    
    if(AUCimp$AUC[j] != 1){
      AUCimp$logit.AUC[j] <- logit(AUCimp$AUC[j])
      AUCimp$logit.se[j]  <- (logit(AUCimp$AUC[j]) - logit(AUCimp$LL[j]))/1.96
      AUCimp$logit.var[j] <- AUCimp$logit.se[j]^2
    } else{
      AUCimp$logit.AUC[j] <- logit(0.999)
      AUCimp$logit.se[j]  <- (logit(0.999) - logit(AUCimp$LL[j]))/1.96
      AUCimp$logit.var[j] <- AUCimp$logit.se[j]^2
    }
    AUCimp$Mal[j] <- nrow(Df[Df$y == 1 & Df$imp == j,])
  }
  
  # Combine results with Rubin's rule
  AUCcombined <- matrix(ncol = 3, nrow = 1)
  colnames(AUCcombined) <- c('logit.AUC', 'logit.LL', 'logit.UL')
  AUCcombined <- data.frame(AUCcombined)
  
  AUCcombined$logit.AUC <- mean(AUCimp$logit.AUC)
  WithinVar <- mean(AUCimp$logit.var)
  BetweenVar <- var(AUCimp$logit.AUC)
  PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
  AUCcombined$PooledSE <- sqrt(PooledVar)
  AUCcombined$logit.LL <- AUCcombined$logit.AUC - 1.96*AUCcombined$PooledSE
  AUCcombined$logit.UL <- AUCcombined$logit.AUC + 1.96*AUCcombined$PooledSE
  
  # Transform back to original scale
  AUCcombined$AUC <- inv.logit(AUCcombined$logit.AUC)
  AUCcombined$LL <- inv.logit(AUCcombined$logit.LL)
  AUCcombined$UL <- inv.logit(AUCcombined$logit.UL)
  
  return(structure(list(Performance = AUCcombined, data = Df, imp = AUCimp)))
  
}

Polcal <- function(outcome,k,p,LP,r=1,estimates=FALSE,dfr=2,plotoverall=TRUE,datapoints=TRUE,smoothing=TRUE,smoothpar=1,intercept=FALSE,slope=FALSE,test=FALSE){
  
  # NOTE: This function is written for a multinomial outcome with three categories.
  #       	  If there are more than three categories, the code can be easily adapted.
  #       	  Comments are added to state where and how the code should be changed if needed.
  
  
  ################################################################################### outcome 	column vector containing the outcome for every case, with values 1 to k (i.c. k=3)
  # k 		number of outcome categories (i.c. 3)
  # p 		matrix with the probabilities of the prediction model, ordered from prob(cat. 1) to 
  # 		prob(cat. k)
  # LP		matrix with all the linear predictors with respect to the chosen reference category, 
  #		ordered (e.g. LP2vs1 and LP3vs1)
  # r 		reference category (default: category 1)
  # estimates 	indicates whether the coefficients of the parametric recalibration framework are desired 
  #		(default=FALSE)
  # dfr 		degrees of freedom for the non-parametric calibration (default=2)
  # plotoverall 	indicates whether overall (non-)parametric calibration plots are constructed 
  #		(default=TRUE)
  # datapoints	indicates whether the individual datapoints are shown in the overall (non-)parametric 
  #		calibration plots (default = TRUE)
  # smoothing 	indicates whether a smoothed line (using cubic splines) is added to the calibration plots 
  #		(default=TRUE)
  # smoothpar 	smoothing parameter for the smoothed line (default=1)
  # intercept 	indicates whether calibration intercepts are desired (default=FALSE)
  # slope 		indicates whether calibration slopes are desired (default=FALSE)
  # test 		indicates whether statistical tests for calibration are desired (default=FALSE)
  ##################################################################################
  
  # for this function you need to use library (VGAM)
  library(VGAM)
  library(bayesm) # necessary for testing
  
  # checks
  if(k != length(table(outcome))){stop('The number of categories in outcome does not equal the specified number of categories.')}
  if(dim(p)[2]!=k){stop('The number of columns in p does not equal the specified number of categories.')}
  if(dim(LP)[2]!=k-1){stop('The number of columns in LP does not equal the specified number of categories minus 1.')}
  if(! r %in% 1:k){stop('The given reference category (r) is not possible.')}     
  if(!is.matrix(p)){stop('p is not a matrix.')}
  if(!is.matrix(LP)){stop('LP is not a matrix.')}
  if(isTRUE(plotoverall) && !isTRUE(datapoints) && !isTRUE(smoothing)){stop('For the overall (non-)parametric calibration plots either datapoints or smoothed lines should be requested.')}
  
  # if tests for perfect calibration are requested, automatically calibration intercepts and calibration slopes 
  # are given
  if(isTRUE(test)){intercept<-slope<-TRUE}
  
  # probabilities
  probs <- split(p,col(p))    
  
  # linear predictors necessary for non-parametric calibration plot - give a name to each linear predictor 
  # seperately
  lps <- split(LP,col(LP))
  for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps[[i]]))}
  
  
  ###############################################
  # parametric logistic recalibration framework 
  # cf. section 2.2.1.                          
  ###############################################
  
  # reference category r
  # LP = matrix with linear predictors
  
  fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
  if(isTRUE(estimates)){est<-coefficients(fitp)
  names(est) <- paste('EST',names(est),sep='.')}
  
  
  #######################################################
  # non-parametric recalibration framework (using df=2) 
  # cf. section 2.2.1.                                  
  #######################################################
  
  # reference category r
  # lpi = ith linear predictor
  # for k outcome categories, there are k-1 linear predictors and the code should be adapted to:
  # fitnp <- vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+...+s(lpk-1,df=dfr),family=multinomial(refLevel=r))
  
  fitnp<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+s(lp3, df=dfr)+s(lp4, df=dfr),family=multinomial(refLevel=r), bf.maxit=100)  
  
  ###############################################                  
  # Separate (non-)parametric calibration plots
  ###############################################
  
  windows()
  par(mfrow=c(ceiling(k/2),2))
  

  
  # non-parametric calibration plot 
  # cf. section 2.2.2.              
  ###################################
  windows()
  par(mfrow=c(ceiling(k/2),2))
  for(i in 1:k){p <- unlist(probs[[i]])
  if(isTRUE(smoothing)){color<-'grey'}else{color<-1+i}
  matplot(p,fitted(fitnp)[,i],type="p",pch=i,col=color,lwd=1,ylab="",xlab="",xlim=0:1,ylim=0:1)
  par(new=T)
  ref <- rbind(c(0,0),c(1,1))
  matplot(ref,ref,type="l",col=1,lwd=2,ylab="Observed proportions",xlab="Predicted probabilities",xlim=0:1,ylim=0:1)
  # smoother for calibration plots 
  ##################################
  # a = smoothing parameter
  if(isTRUE(smoothing)){
    a = smoothpar
    points(smooth.spline(p, fitted(fitnp)[,i],spar=a), type="l", col=(1+i), lwd = 4)}
  # legend
  legende <- c(paste("cat ", i, sep = ""))
  legend(x=0.6, y=(0.2),col=(1+i),lty =1,legend=legende)
  title(main = "Non-parametric calibration plot")
  par(new=F)}
  
  
  #############################################            
  # Overall (non-)parametric calibration plot 
  #############################################
  
  if(isTRUE(plotoverall)){
    windows()
    

    
    # non-parametric calibration plot 
    # cf. section 2.2.2.              
    ###################################
    
    windows()
    if(isTRUE(datapoints)){for(i in 1:k){p <- unlist(probs[[i]])
    matplot(p,fitted(fitnp)[,i],type="p",pch=i,col=(1+i),lwd=1,ylab="",xlab="",xlim=0:1,ylim=0:1)
    par(new=T)}}
    ref <- rbind(c(0,0),c(1,1))
    matplot(ref,ref,type="l",col=1,lwd=2,ylab="Observed proportions",xlab="Predicted  probabilities",xlim=0:1,ylim=0:1)
    # smoother for calibration plots 
    ##################################
    # a = smoothing parameter
    if(isTRUE(smoothing)){a = smoothpar
    for(i in 1:k){p <- unlist(probs[[i]])
    points(smooth.spline(p, fitted(fitnp)[,i],df=8), type="l", col=cols[i], lwd = 4)}}
    # legend
    for(i in 1:k){if(i <= 2){legende <- c("cat 1","cat 2")}
      if(i > 2){legende <- c("benign", "borderline", "stage I", "stage II-IV", "secondary metastatic")}}
    legend(x=0.6, y=(0.20+(k-3)*0.05),col=c("#440154FF", "#39568CFF","#1F968BFF","#73D055FF","#FDE725FF"),lty =1,legend=legende, lwd=4)
    title(main = "mixed-effects LR with colscore: flexible multinomial calibration plot", cex.main=1)
    par(new=F)}
  
  
  ########################################
  # estimation of calibration intercepts 
  # cf. section 2.2.3. and 2.2.4.        
  ########################################
  
  if(isTRUE(intercept)){int<-vgam(outcome~1,offset=LP,family=multinomial(refLevel=r))
  coeffint<-coefficients(int)
  se<-sqrt(diag(vcov(int)))
  ci1i <- cbind(LL1 = coeffint[1] - qnorm(0.975) * se[1], UL1 = coeffint[1] + qnorm(0.975) * se[1])
  ci2i <- cbind(LL2 = coeffint[2] - qnorm(0.975) * se[2], UL2 = coeffint[2] + qnorm(0.975) * se[2])
  ci3i <- cbind(LL1 = coeffint[3] - qnorm(0.975) * se[3], UL1 = coeffint[3] + qnorm(0.975) * se[3])
  ci4i <- cbind(LL1 = coeffint[4] - qnorm(0.975) * se[4], UL1 = coeffint[4] + qnorm(0.975) * se[4])
  
  estint <- c(coeffint[1],ci1i,coeffint[2],ci2i, coeffint[3], ci3i, coeffint[4], ci4i)
  names(estint) <- paste('CALINT',c("int1","LLint1","ULint1","int2","LLint2","ULint2", "int3", "LLint3", "ULint3", "int4", "LLint4", "ULint4"),sep='.')
  }
  
  
  ####################################
  # estimation of calibration slopes 
  # cf. section 2.2.3. and 2.2.4.    
  ####################################
  
  # we used constraints to fix some coefficients to zero as appropriate
  # for k outcome categories this code should be changed to:
  # i <- diag(k-1)
  # i2 <- cbind(c(1,rep(0,k-2)))
  # i3 <- cbind(c(0,1,rep(0,k-1)))
  # i4 <- cbind(c(0,0,1,rep(0,k-2)))
  # ... (ij <- cbind(c(rep(0,j-2),1,rep(0,k-j)))
  # ik <- cbind(c(rep(0,k-2),1))
  # clist<-list("(Intercept)"=i,"lp1"=i2,"lp2"=i3,...,"lpk-1"=ik)
  # slopes<-vgam(outcome~lp1+lp2+...+lpk-1,family=multinomial(refLevel=r),constraints=clist)
  
  if(isTRUE(slope)){
    i<-diag(k-1)
    i2<-cbind(c(1, rep(0, k-2)))
    i3<-cbind(c(0,1, rep(0, k-3)))
    i4 <- cbind(c(0,0,1,0))
    i5 <- cbind(c(0,0,0,1))
    clist<-list("(Intercept)"=i,"lp1"=i2,"lp2"=i3, "lp3"=i4, "lp4"=i5)
    slopes<-vgam(outcome~lp1+lp2+lp3+lp4,family=multinomial(refLevel=r),constraints=clist)
    coeffslopes<-coefficients(slopes)[k:length(coefficients(slopes))]
    se<-sqrt(diag(vcov(slopes)))
    ci1s <- cbind(LL1 = coeffslopes[1] - qnorm(0.975) * se[5], UL1 = coeffslopes[1] + qnorm(0.975) * se[5])
    ci2s <- cbind(LL2 = coeffslopes[2] - qnorm(0.975) * se[6], UL2 = coeffslopes[2] + qnorm(0.975) * se[6])
    ci3s <- cbind(LL1 = coeffslopes[3] - qnorm(0.975) * se[7], UL1 = coeffslopes[3] + qnorm(0.975) * se[7])
    ci4s <- cbind(LL1 = coeffslopes[4] - qnorm(0.975) * se[8], UL1 = coeffslopes[4] + qnorm(0.975) * se[8])
    
    estslopes <- c(coeffslopes[1],ci1s,coeffslopes[2],ci2s, coeffslopes[3], ci3s, coeffslopes[4], ci4s)
    names(estslopes) <- paste('CALSLOPES',c('lp1','LLlp1','ULlp1','lp2','LLlp2','ULlp2','lp3','LLlp3','ULlp3','lp4','LLlp4','ULlp4'),sep='.')}
  
  
  #################################
  # calibration testing          
  # cf. section 2.2.3. and 2.2.4. 
  #################################
  
  # this code requires the bayesm library developed by Peter Rossi
  
  if(isTRUE(test)){
    
    # -2 log-likelihood of model without adaptations
    # for k outcome categories this code should be changed to:
    # alphas <- rep(0,k-1) #(i.e. all intercepts zero)
    # beta1 <- c(1,rep(0,k-2)) #(i.e. first linear predictor for first equation)
    # beta2 <- c(0,1,rep(0,k-3)) #(i.e. second linear predictor for second equation)      
    # betaj <- c(rep(0,j-1),1,rep(0,k-1-j)) #(i.e. jth linear predictor for jth equation)
    # betak <- c(rep(0,k-2),1) #(i.e. kth linear predictor for kth equation)
    # parametersk <- c(alphas, beta1, beta2, ..., betak)
    
    parametersk <- c(0,0,1,0,0,1) #c(alpha1,alpha2,b22,b23,b32,b33)
    Xdk=LP
    x <- createX(p=k,na=0,nd=k-1,Xa=NULL,Xd=Xdk,INT=TRUE,DIFF=FALSE,base=1)
    deviancewithout <- -2*llmnl(parametersk,outcome,x)
    names(deviancewithout)<-c('original deviance')
    
    devint <- deviance(int)
    names(devint)<-c('intercept deviance')
    devslopes <- deviance(slopes)
    names(devslopes)<-c('slopes deviance')
    
    # overall calibration (i.e. calibration intercepts and slopes) 
    ################################################################
    
    poverall<- pchisq(deviancewithout - devslopes, df = 2*(k-1), lower.tail = FALSE)
    
    # calibration intercepts 
    ##########################
    
    pint<- pchisq(deviancewithout - devint, df = k-1, lower.tail = FALSE)
    
    # calibration slopes 
    ######################
    
    pslopes<- pchisq(devint - devslopes, df = k-1, lower.tail = FALSE)
    names(poverall)<-c('p overall')
    names(pint)<-c('p int')
    names(pslopes)<-c('p slopes')}
  
  # Printing of results
  # The probabilities of calibration intercepts and slopes are only shown when the hypothesis of perfect 
  # calibration is rejected.
  
  results<-list(if(isTRUE(estimates)){est}else{'Not requested'},if(isTRUE(intercept)){estint}else{'Not requested'},if(isTRUE(slope)){estslopes}else{'Not requested'},if(isTRUE(test)){c(deviancewithout,devint,devslopes)}else{'Not requested'},if(isTRUE(test)){c(poverall,if(poverall<0.05){c(pint,pslopes)})}else{'Not requested'})
  names(results)<-c("Coefficients of parametric recalibration framework","Calibration Intercepts with 95% CI","Calibration Slopes with 95% CI","Deviances","P-values")
  n <- 1:5
  selection <- c(isTRUE(estimates),isTRUE(intercept),isTRUE(slope),isTRUE(test),isTRUE(test))
  results[n[selection]]
  
}



#clinical utility
DataWinBugs.imp <- function(pred, outcome, center, data, imp,
                            sequence = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  NRimp <- length(unique(Df$imp))
  
  centers <- unique(Df$center)
  
  ConfusionList <- list()
  for(i in 1:length(sequence)){
    ConfusionList[[i]] <- list()
  }
  
  ConfusionImp <- list()
  for(i in 1:NRimp){
    ConfusionImp[[i]] <- list()
  }
  
  for(k in 1:NRimp){
    cat("Imputation", k, "of", NRimp)
    for(i in 1:length(sequence)){
      threshold <- sequence[i]
      
      Confusion <- matrix(nrow = length(centers), ncol = 8)
      Confusion <- data.frame(Confusion)
      colnames(Confusion) <- c('Center', 'CutOff', 'TN', 'TP', 'FP', 'FN', 'cases', 'controls')
      Confusion$CutOff <- threshold
      
      for(j in seq_along(centers)){
        
        Confusion$Center[j] <- centers[j]
        
        CM <- confusion.matrix(obs = Df$y[Df$center == centers[j] & Df$imp == k], pred = Df$p[Df$center == centers[j] & Df$imp == k], threshold = threshold)
        Confusion$TN[j] <- CM[1,1]
        Confusion$TP[j] <- CM[2,2]
        Confusion$FP[j] <- CM[2,1]
        Confusion$FN[j] <- CM[1,2]
        
        Confusion$cases <- Confusion$TP + Confusion$FN
        Confusion$controls <- Confusion$TN + Confusion$FP
        
        Confusion$n <- Confusion$cases + Confusion$controls
        Confusion$NB <- Confusion$TP / Confusion$n - Confusion$FP / Confusion$n * (threshold / (1 - threshold))
        
      }
      
      ConfusionList[[i]] <- Confusion
    }
    cat("rbindlist toepassen")
    ConfusionImp[[k]] <- rbindlist(ConfusionList, fill = TRUE)
  }
  
  ConfusionLong <- rbindlist(ConfusionImp, fill = TRUE)
  ConfusionSum <- summaryBy(cbind(TP, TN, cases, controls, FP, FN, NB) ~ cbind(CutOff, Center), data = ConfusionLong, FUN = mean)
  
  return(structure(list(Results = ConfusionSum)))
}

## Without multiple imputation
DataWinBugs <- function(pred, outcome, center, data, 
                        sequence = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  ConfusionList <- list()
  
  for(i in 1:length(sequence)){
    ConfusionList[[i]] <- list()
  }
  
  
  for(i in 1:length(sequence)){
    threshold <- sequence[i]
    
    Confusion <- matrix(nrow = length(centers), ncol = 8)
    Confusion <- data.frame(Confusion)
    colnames(Confusion) <- c('Center', 'CutOff', 'TN', 'TP', 'FP', 'FN', 'cases', 'controls')
    Confusion$CutOff <- threshold
    
    for(j in seq_along(centers)){
      
      Confusion$Center[j] <- centers[j]
      
      CM <- confusion.matrix(obs = Df$y[Df$center == centers[j]], pred = Df$p[Df$center == centers[j]], threshold = threshold)
      Confusion$TN[j] <- CM[1,1]
      Confusion$TP[j] <- CM[2,2]
      Confusion$FP[j] <- CM[2,1]
      Confusion$FN[j] <- CM[1,2]
      
      Confusion$cases <- Confusion$TP + Confusion$FN
      Confusion$controls <- Confusion$TN + Confusion$FP
      
      Confusion$n <- Confusion$cases + Confusion$controls
      Confusion$NB <- Confusion$TP / Confusion$n - Confusion$FP / Confusion$n * (threshold / (1 - threshold))
      
    }
    
    ConfusionList[[i]] <- Confusion
    
    
  }
  return(structure(list(Results = ConfusionList)))
}

#Sensitivity/Specificity analyses

## With multiple imputation: fixed threshold
SS.imp <- function(pred, outcome, threshold, center, imp, data, method.MA = "REML"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  threshold = threshold
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  
  Sensimp <- list()
  for(i in 1:NRimp){
    Sensimp[[i]] <- list()
  }
  
  
  # Sensitivity and specificity per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    Senscenter <- matrix(ncol = 7, nrow = length(centers))
    colnames(Senscenter) <- c("Center", "TP", "TN", "FP", "FN", "Sensitivity", "Specificity")
    Senscenter <- data.frame(Senscenter)
    Senscenter$Center <- centers
    
    for(i in seq_along(centers)){
      predicted_values <- ifelse(Df$p[Df$center == centers[i] & Df$imp == j] >= threshold, 1, 0)
      actual_values <- Df$y[Df$center == centers[i] & Df$imp == j]
      conf_matrix <- table(factor(predicted_values, levels = c(0, 1)), factor(actual_values, levels = c(0, 1)))
      Senscenter$TN[i]     <- conf_matrix[1,1]
      Senscenter$TP[i]     <- conf_matrix[2,2]
      Senscenter$FN[i]     <- conf_matrix[1,2]
      Senscenter$FP[i]     <- conf_matrix[2,1]
      
      if(any(Senscenter[i, 2:5] == 0)){
        Senscenter$TN[i]     <- conf_matrix[1,1] + 0.1
        Senscenter$TP[i]     <- conf_matrix[2,2] + 0.1
        Senscenter$FN[i]     <- conf_matrix[1,2] + 0.1
        Senscenter$FP[i]     <- conf_matrix[2,1] + 0.1
      }
      
    }
    
    Senscenter$Sensitivity  <- Senscenter$TP / (Senscenter$TP + Senscenter$FN)
    Senscenter$se_sens      <- sqrt((Senscenter$Sensitivity * (1 - Senscenter$Sensitivity))/(Senscenter$TP + Senscenter$FN))
    Senscenter$Specificity  <- Senscenter$TN / (Senscenter$TN + Senscenter$FP)
    Senscenter$se_spec      <- sqrt((Senscenter$Specificity * (1 - Senscenter$Specificity))/(Senscenter$TN + Senscenter$FP))
    
    # Logit transformation
    Senscenter$logit.sens     <- logit(Senscenter$Sensitivity)
    Senscenter$logit.se.sens  <- sqrt(1/Senscenter$TP + 1/Senscenter$FN)
    Senscenter$logit.var.sens <- Senscenter$logit.se.sens^2
    
    Senscenter$logit.spec     <- logit(Senscenter$Specificity)
    Senscenter$logit.se.spec  <- sqrt(1/Senscenter$TN + 1/Senscenter$FP)
    Senscenter$logit.var.spec <- Senscenter$logit.se.spec^2
    
    Sensimp[[j]] <- Senscenter
    
  }
  
  SensimpLong <- rbindlist(Sensimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  Senscombined <- matrix(ncol = 7, nrow = length(centers))
  colnames(Senscombined) <- c('Center', 'logit.sens', 'logit.LL.sens', 'logit.UL.sens', 'logit.spec', 'logit.LL.spec', 'logit.UL.spec')
  Senscombined <- data.frame(Senscombined)
  Senscombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    
    Senscombined$TP[i] <- mean(SensimpLong$TP[SensimpLong$Center == centers[i]])
    Senscombined$TN[i] <- mean(SensimpLong$TN[SensimpLong$Center == centers[i]])
    Senscombined$FP[i] <- mean(SensimpLong$FP[SensimpLong$Center == centers[i]])
    Senscombined$FN[i] <- mean(SensimpLong$FN[SensimpLong$Center == centers[i]])
    
    # Sensitivity
    Senscombined$logit.sens[i] <- mean(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    WithinVar  <- mean(SensimpLong$logit.var.sens[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    PooledVar  <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.sens[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.sens[i] <- Senscombined$logit.sens[i] - 1.96*Senscombined$PooledSE.sens[i]
    Senscombined$logit.UL.sens[i] <- Senscombined$logit.sens[i] + 1.96*Senscombined$PooledSE.sens[i]
    
    # Specificity
    Senscombined$logit.spec[i] <- mean(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    WithinVar <- mean(SensimpLong$logit.var.spec[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.spec[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.spec[i] <- Senscombined$logit.spec[i] - 1.96*Senscombined$PooledSE.spec[i]
    Senscombined$logit.UL.spec[i] <- Senscombined$logit.spec[i] + 1.96*Senscombined$PooledSE.spec[i]
  }
  
  # Transform back to original scale
  Senscombined$Sensitivity <- inv.logit(Senscombined$logit.sens)
  Senscombined$LL.sens     <- inv.logit(Senscombined$logit.LL.sens)
  Senscombined$UL.sens     <- inv.logit(Senscombined$logit.UL.sens)
  
  Senscombined$Specificity <- inv.logit(Senscombined$logit.spec)
  Senscombined$LL.spec     <- inv.logit(Senscombined$logit.LL.spec)
  Senscombined$UL.spec     <- inv.logit(Senscombined$logit.UL.spec)
  
  
  Sensoverall <- matrix(nrow = 1, ncol = 7)
  colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
  Sensoverall <- data.frame(Sensoverall)
  Sensoverall$Center <- "Overall"
  
  
  ## Bivariate random-effects model
  p1 <- 2*length(centers)
  p2 <- 4*length(centers)
  Sensbi <- Senscombined[, c("Center", "logit.sens", "logit.spec", "PooledSE.sens", "PooledSE.spec")]
  Sensbi.long <- melt(Sensbi, id.vars = "Center")
  Sensbi2 <- cbind(Sensbi.long[1:p1,], Sensbi.long[(p1 + 1):p2,])
  colnames(Sensbi2) <- c("Center", "SensSpec", "Value", "Centre", "Pooled", "SE")
  Sensbi2$VAR <- Sensbi2$SE^2
  Sensbi2$SensSpec <- factor(Sensbi2$SensSpec)
  fit.RE = rma.mv(yi = Value, V = VAR, mods = ~ SensSpec-1, random = ~ SensSpec-1 | Center, struct="UN", data=Sensbi2)
  Sensoverall$Sens <- inv.logit(coef(fit.RE)[1])
  Sensoverall$LL.sens <- inv.logit(fit.RE$ci.lb[1])
  Sensoverall$UL.sens <- inv.logit(fit.RE$ci.ub[1])
  Sensoverall$Spec <- inv.logit(coef(fit.RE)[2])
  Sensoverall$LL.spec <- inv.logit(fit.RE$ci.lb[2])
  Sensoverall$UL.spec <- inv.logit(fit.RE$ci.ub[2])
  
  Sensoverall$Threshold   <- threshold
  Sensoverall$Sensitivity <- paste0(format(round(Sensoverall$Sens, 3), nsmall = 3), " (", format(round(Sensoverall$LL.sens, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.sens, 3), nsmall = 3), ")")
  Sensoverall$Specificity <- paste0(format(round(Sensoverall$Spec, 3), nsmall = 3), " (", format(round(Sensoverall$LL.spec, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.spec, 3), nsmall = 3), ")")
  
  return(structure(list(OverallPer = Sensoverall[, 8:10], CenterPer = Senscombined, IncludedCenter = centers, data = Df)))
}

## With multiple imputation: fixed sensitivity
FixedSens.imp <- function(pred, outcome, Sensitivity, center, imp, data, method.MA = "REML"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  Sensitivity = Sensitivity
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  
  Sensimp <- list()
  for(i in 1:NRimp){
    Sensimp[[i]] <- list()
  }
  
  
  # Sensitivity and specificity per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    Senscenter <- matrix(ncol = 3, nrow = length(centers))
    colnames(Senscenter) <- c("Center", "Sensitivity", "Specificity")
    Senscenter <- data.frame(Senscenter)
    Senscenter$Center <- centers
    Senscenter$Sensitivity <- Sensitivity
    
    for(i in seq_along(centers)){
      threshold <- quantile(Df$p[Df$center == centers[i] & Df$imp == j & Df$y == 1], probs = 1 - Sensitivity)
      
      predicted_values <- ifelse(Df$p[Df$center == centers[i] & Df$imp == j] >= threshold, 1, 0)
      actual_values <- Df$y[Df$center == centers[i] & Df$imp == j]
      conf_matrix <- table(factor(predicted_values, levels = c(0, 1)), factor(actual_values, levels = c(0, 1)))
      Senscenter$TN[i]     <- conf_matrix[1,1]
      Senscenter$TP[i]     <- conf_matrix[2,2]
      Senscenter$FN[i]     <- conf_matrix[1,2]
      Senscenter$FP[i]     <- conf_matrix[2,1]
      
      if(any(Senscenter[i, 4:7] == 0)){
        Senscenter$TN[i]     <- conf_matrix[1,1] + 0.1
        Senscenter$TP[i]     <- conf_matrix[2,2] + 0.1
        Senscenter$FN[i]     <- conf_matrix[1,2] + 0.1
        Senscenter$FP[i]     <- conf_matrix[2,1] + 0.1
      }
      
    }
    
    Senscenter$Sensitivity  <- Senscenter$TP / (Senscenter$TP + Senscenter$FN)
    Senscenter$se_sens      <- sqrt((Senscenter$Sensitivity * (1 - Senscenter$Sensitivity))/(Senscenter$TP + Senscenter$FN))
    Senscenter$Specificity  <- Senscenter$TN / (Senscenter$TN + Senscenter$FP)
    Senscenter$se_spec      <- sqrt((Senscenter$Specificity * (1 - Senscenter$Specificity))/(Senscenter$TN + Senscenter$FP))
    
    # Logit transformation
    Senscenter$logit.sens     <- logit(Senscenter$Sensitivity)
    Senscenter$logit.se.sens  <- sqrt(1/Senscenter$TP + 1/Senscenter$FN)
    Senscenter$logit.var.sens <- Senscenter$logit.se.sens^2
    
    Senscenter$logit.spec     <- logit(Senscenter$Specificity)
    Senscenter$logit.se.spec  <- sqrt(1/Senscenter$TN + 1/Senscenter$FP)
    Senscenter$logit.var.spec <- Senscenter$logit.se.spec^2
    
    Sensimp[[j]] <- Senscenter
    
  }
  
  SensimpLong <- rbindlist(Sensimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  Senscombined <- matrix(ncol = 7, nrow = length(centers))
  colnames(Senscombined) <- c('Center', 'logit.sens', 'logit.LL.sens', 'logit.UL.sens', 'logit.spec', 'logit.LL.spec', 'logit.UL.spec')
  Senscombined <- data.frame(Senscombined)
  Senscombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    
    Senscombined$TP[i] <- mean(SensimpLong$TP[SensimpLong$Center == centers[i]])
    Senscombined$TN[i] <- mean(SensimpLong$TN[SensimpLong$Center == centers[i]])
    Senscombined$FP[i] <- mean(SensimpLong$FP[SensimpLong$Center == centers[i]])
    Senscombined$FN[i] <- mean(SensimpLong$FN[SensimpLong$Center == centers[i]])
    
    # Sensitivity
    Senscombined$logit.sens[i] <- mean(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    WithinVar  <- mean(SensimpLong$logit.var.sens[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    PooledVar  <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.sens[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.sens[i] <- Senscombined$logit.sens[i] - 1.96*Senscombined$PooledSE.sens[i]
    Senscombined$logit.UL.sens[i] <- Senscombined$logit.sens[i] + 1.96*Senscombined$PooledSE.sens[i]
    
    # Specificity
    Senscombined$logit.spec[i] <- mean(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    WithinVar <- mean(SensimpLong$logit.var.spec[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.spec[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.spec[i] <- Senscombined$logit.spec[i] - 1.96*Senscombined$PooledSE.spec[i]
    Senscombined$logit.UL.spec[i] <- Senscombined$logit.spec[i] + 1.96*Senscombined$PooledSE.spec[i]
  }
  
  # Transform back to original scale
  Senscombined$Sensitivity <- inv.logit(Senscombined$logit.sens)
  Senscombined$LL.sens     <- inv.logit(Senscombined$logit.LL.sens)
  Senscombined$UL.sens     <- inv.logit(Senscombined$logit.UL.sens)
  
  Senscombined$Specificity <- inv.logit(Senscombined$logit.spec)
  Senscombined$LL.spec     <- inv.logit(Senscombined$logit.LL.spec)
  Senscombined$UL.spec     <- inv.logit(Senscombined$logit.UL.spec)
  
  
  Sensoverall <- matrix(nrow = 1, ncol = 7)
  colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
  Sensoverall <- data.frame(Sensoverall)
  Sensoverall$Center <- "Overall"
  
  
  ## Bivariate random-effects model
  p1 <- 2*length(centers)
  p2 <- 4*length(centers)
  Sensbi <- Senscombined[, c("Center", "logit.sens", "logit.spec", "PooledSE.sens", "PooledSE.spec")]
  Sensbi.long <- melt(Sensbi, id.vars = "Center")
  Sensbi2 <- cbind(Sensbi.long[1:p1,], Sensbi.long[(p1 + 1):p2,])
  colnames(Sensbi2) <- c("Center", "SensSpec", "Value", "Centre", "Pooled", "SE")
  Sensbi2$VAR <- Sensbi2$SE^2
  Sensbi2$SensSpec <- factor(Sensbi2$SensSpec)
  fit.RE = rma.mv(yi = Value, V = VAR, mods = ~ SensSpec-1, random = ~ SensSpec | Center, struct="UN", data=Sensbi2)
  Sensoverall$Sens <- inv.logit(coef(fit.RE)[1])
  Sensoverall$LL.sens <- inv.logit(fit.RE$ci.lb[1])
  Sensoverall$UL.sens <- inv.logit(fit.RE$ci.ub[1])
  Sensoverall$Spec <- inv.logit(coef(fit.RE)[2])
  Sensoverall$LL.spec <- inv.logit(fit.RE$ci.lb[2])
  Sensoverall$UL.spec <- inv.logit(fit.RE$ci.ub[2])
  
  Sensoverall$Threshold   <- threshold
  Sensoverall$Sensitivity <- paste0(format(round(Sensoverall$Sens, 3), nsmall = 3), " (", format(round(Sensoverall$LL.sens, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.sens, 3), nsmall = 3), ")")
  Sensoverall$Specificity <- paste0(format(round(Sensoverall$Spec, 3), nsmall = 3), " (", format(round(Sensoverall$LL.spec, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.spec, 3), nsmall = 3), ")")
  
  return(structure(list(OverallPer = Sensoverall[, 8:10], CenterPer = Senscombined, IncludedCenter = centers, data = Df)))
}

## With multiple imputation: fixed specificity
FixedSpec.imp <- function(pred, outcome, Specificity, center, imp, data, method.MA = "REML"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  Specificity = Specificity
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  
  Sensimp <- list()
  for(i in 1:NRimp){
    Sensimp[[i]] <- list()
  }
  
  
  # Sensitivity and specificity per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    Senscenter <- matrix(ncol = 3, nrow = length(centers))
    colnames(Senscenter) <- c("Center", "Sensitivity", "Specificity")
    Senscenter <- data.frame(Senscenter)
    Senscenter$Center <- centers
    Senscenter$Specificity <- Specificity
    
    for(i in seq_along(centers)){
      threshold <- quantile(Df$p[Df$center == centers[i] & Df$imp == j & Df$y == 0], probs = Specificity)
      
      predicted_values <- ifelse(Df$p[Df$center == centers[i] & Df$imp == j] >= threshold, 1, 0)
      actual_values <- Df$y[Df$center == centers[i] & Df$imp == j]
      conf_matrix <- table(factor(predicted_values, levels = c(0, 1)), factor(actual_values, levels = c(0, 1)))
      Senscenter$TN[i]     <- conf_matrix[1,1]
      Senscenter$TP[i]     <- conf_matrix[2,2]
      Senscenter$FN[i]     <- conf_matrix[1,2]
      Senscenter$FP[i]     <- conf_matrix[2,1]
      
      if(any(Senscenter[i, 4:7] == 0)){
        Senscenter$TN[i]     <- conf_matrix[1,1] + 0.1
        Senscenter$TP[i]     <- conf_matrix[2,2] + 0.1
        Senscenter$FN[i]     <- conf_matrix[1,2] + 0.1
        Senscenter$FP[i]     <- conf_matrix[2,1] + 0.1
      }
      
    }
    
    Senscenter$Sensitivity  <- Senscenter$TP / (Senscenter$TP + Senscenter$FN)
    Senscenter$se_sens      <- sqrt((Senscenter$Sensitivity * (1 - Senscenter$Sensitivity))/(Senscenter$TP + Senscenter$FN))
    Senscenter$Specificity  <- Senscenter$TN / (Senscenter$TN + Senscenter$FP)
    Senscenter$se_spec      <- sqrt((Senscenter$Specificity * (1 - Senscenter$Specificity))/(Senscenter$TN + Senscenter$FP))
    
    # Logit transformation
    Senscenter$logit.sens     <- logit(Senscenter$Sensitivity)
    Senscenter$logit.se.sens  <- sqrt(1/Senscenter$TP + 1/Senscenter$FN)
    Senscenter$logit.var.sens <- Senscenter$logit.se.sens^2
    
    Senscenter$logit.spec     <- logit(Senscenter$Specificity)
    Senscenter$logit.se.spec  <- sqrt(1/Senscenter$TN + 1/Senscenter$FP)
    Senscenter$logit.var.spec <- Senscenter$logit.se.spec^2
    
    Sensimp[[j]] <- Senscenter
    
  }
  
  SensimpLong <- rbindlist(Sensimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  Senscombined <- matrix(ncol = 7, nrow = length(centers))
  colnames(Senscombined) <- c('Center', 'logit.sens', 'logit.LL.sens', 'logit.UL.sens', 'logit.spec', 'logit.LL.spec', 'logit.UL.spec')
  Senscombined <- data.frame(Senscombined)
  Senscombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    
    Senscombined$TP[i] <- mean(SensimpLong$TP[SensimpLong$Center == centers[i]])
    Senscombined$TN[i] <- mean(SensimpLong$TN[SensimpLong$Center == centers[i]])
    Senscombined$FP[i] <- mean(SensimpLong$FP[SensimpLong$Center == centers[i]])
    Senscombined$FN[i] <- mean(SensimpLong$FN[SensimpLong$Center == centers[i]])
    
    # Sensitivity
    Senscombined$logit.sens[i] <- mean(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    WithinVar  <- mean(SensimpLong$logit.var.sens[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.sens[SensimpLong$Center == centers[i]])
    PooledVar  <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.sens[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.sens[i] <- Senscombined$logit.sens[i] - 1.96*Senscombined$PooledSE.sens[i]
    Senscombined$logit.UL.sens[i] <- Senscombined$logit.sens[i] + 1.96*Senscombined$PooledSE.sens[i]
    
    # Specificity
    Senscombined$logit.spec[i] <- mean(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    WithinVar <- mean(SensimpLong$logit.var.spec[SensimpLong$Center == centers[i]])
    BetweenVar <- var(SensimpLong$logit.spec[SensimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    Senscombined$PooledSE.spec[i] <- sqrt(PooledVar)
    Senscombined$logit.LL.spec[i] <- Senscombined$logit.spec[i] - 1.96*Senscombined$PooledSE.spec[i]
    Senscombined$logit.UL.spec[i] <- Senscombined$logit.spec[i] + 1.96*Senscombined$PooledSE.spec[i]
  }
  
  # Transform back to original scale
  Senscombined$Sensitivity <- inv.logit(Senscombined$logit.sens)
  Senscombined$LL.sens     <- inv.logit(Senscombined$logit.LL.sens)
  Senscombined$UL.sens     <- inv.logit(Senscombined$logit.UL.sens)
  
  Senscombined$Specificity <- inv.logit(Senscombined$logit.spec)
  Senscombined$LL.spec     <- inv.logit(Senscombined$logit.LL.spec)
  Senscombined$UL.spec     <- inv.logit(Senscombined$logit.UL.spec)
  
  
  Sensoverall <- matrix(nrow = 1, ncol = 7)
  colnames(Sensoverall) <- c('Center', 'Sens', 'LL.sens', 'UL.sens', 'Spec', 'LL.spec', 'UL.spec')
  Sensoverall <- data.frame(Sensoverall)
  Sensoverall$Center <- "Overall"
  
  
  ## Bivariate random-effects model
  p1 <- 2*length(centers)
  p2 <- 4*length(centers)
  Sensbi <- Senscombined[, c("Center", "logit.sens", "logit.spec", "PooledSE.sens", "PooledSE.spec")]
  Sensbi.long <- melt(Sensbi, id.vars = "Center")
  Sensbi2 <- cbind(Sensbi.long[1:p1,], Sensbi.long[(p1 + 1):p2,])
  colnames(Sensbi2) <- c("Center", "SensSpec", "Value", "Centre", "Pooled", "SE")
  Sensbi2$VAR <- Sensbi2$SE^2
  Sensbi2$SensSpec <- factor(Sensbi2$SensSpec)
  fit.RE = rma.mv(yi = Value, V = VAR, mods = ~ SensSpec-1, random = ~ SensSpec | Center, struct="UN", data=Sensbi2)
  Sensoverall$Sens <- inv.logit(coef(fit.RE)[1])
  Sensoverall$LL.sens <- inv.logit(fit.RE$ci.lb[1])
  Sensoverall$UL.sens <- inv.logit(fit.RE$ci.ub[1])
  Sensoverall$Spec <- inv.logit(coef(fit.RE)[2])
  Sensoverall$LL.spec <- inv.logit(fit.RE$ci.lb[2])
  Sensoverall$UL.spec <- inv.logit(fit.RE$ci.ub[2])
  
  Sensoverall$Threshold   <- threshold
  Sensoverall$Sensitivity <- paste0(format(round(Sensoverall$Sens, 3), nsmall = 3), " (", format(round(Sensoverall$LL.sens, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.sens, 3), nsmall = 3), ")")
  Sensoverall$Specificity <- paste0(format(round(Sensoverall$Spec, 3), nsmall = 3), " (", format(round(Sensoverall$LL.spec, 3), nsmall = 3), " - ", format(round(Sensoverall$UL.spec, 3), nsmall = 3), ")")
  
  return(structure(list(OverallPer = Sensoverall[, 8:10], CenterPer = Senscombined[, c(1, 10:15)], IncludedCenter = centers, data = Df)))
}

AUC.IOTA <- function(pred, outcome, center, data, method.MA = "REML", titleGraph = "AUC per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  AUCcenter <- matrix(ncol = 6, nrow = length(centers))
  colnames(AUCcenter) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCcenter <- data.frame(AUCcenter)
  AUCcenter$Center <- centers
  
  # AUC per center
  for(i in seq_along(centers)){
    AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 0], method = "pepe")
    AUCcenter[i, 2]   <- nrow(Df[Df$center == centers[i],])
    AUCcenter[i, 3]   <- round(nrow(Df[Df$y == 1 & Df$center == centers[i],])/nrow(Df[Df$center == centers[i],])*100)
    
    ## Additional part for AUCs of 1
    if(AUCcenter[i, 4] == 1){
      AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 0], method = "newcombe") # Newcombe ipv pepe
    } else{
      AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 0], method = "pepe")
    }
    
    if(AUCcenter$AUC[i] != 1){
      AUCcenter$logit.AUC[i] <- logit(AUCcenter$AUC[i])
      AUCcenter$logit.se[i]  <- (logit(AUCcenter$AUC[i]) - logit(AUCcenter$LL[i]))/1.96
      AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
    } else{
      AUCcenter$logit.AUC[i] <- logit(0.999)
      AUCcenter$logit.se[i]  <- (logit(0.999) - logit(AUCcenter$LL[i]))/1.96
      AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
    }
  }
  
  
  AUCcenter$logit.AUC <- logit(AUCcenter$AUC)
  AUCcenter$logit.se  <- (logit(AUCcenter$AUC) - logit(AUCcenter$LL))/1.96
  AUCcenter <- AUCcenter[order(AUCcenter$SampleSize, decreasing = TRUE),]
  
  AUCoverall <- matrix(nrow = 2, ncol = 6)
  colnames(AUCoverall) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCoverall <- data.frame(AUCoverall)
  AUCoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  AUCoverall$SampleSize <- nrow(Df)
  AUCoverall$Prevalence <- round(nrow(Df[Df$y == 1,])/nrow(Df)*100)
  
  # Meta-analyse voor overall estimate
  fit.RE = rma.uni(AUCcenter$logit.AUC, sei = AUCcenter$logit.se, method = method.MA)
  PI = predict(fit.RE, transf = transf.ilogit)
  
  AUCoverall$AUC[1] <- inv.logit(coef(fit.RE))
  AUCoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  AUCoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  AUCoverall$AUC[2] <- PI$pred
  AUCoverall$LL[2] <- PI$cr.lb
  AUCoverall$UL[2] <- PI$cr.ub
  AUCoverall
  
  NAforest <- matrix(nrow = 1, ncol = 6)
  colnames(NAforest) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  AUC <- rbind(AUCcenter[, 1:6], NAforest, AUCoverall)
  
  # Layout forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(AUC)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- AUC$Center[i]
  }
  RRauc <- 1:nrobs
  for(i in 1:nrobs){
    RRauc[i] <- paste(format(round(AUC$AUC[i], 2), nsmall = 2), " (", format(round(AUC$LL[i], 2), nsmall = 2), " to ", format(round(AUC$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- AUC$SampleSize[i]
  }
  
  Labels <- c('Centre', 'AUC (95% CI)', 'N')
  Combined <- cbind(RRcenter, RRauc, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  AUCna <- rbind(NAforest, NAforest, AUC)
  
  return(structure(list(Performance = AUCoverall, ModelFit = fit.RE, AUCcenters = AUCcenter, IncludedCenters = centers, dataPlot = AUCna, Plot = Tabletext, data = Df)))
  
}

#RE.ValProb2
RE.ValProb2 <- function (p, y, center, data, CalibrLines = c("overall", "centers", 
                                                             "both"), LogCal = T, flexible = F, dostats = T, statloc = c(0, 
                                                                                                                         0.85), legendloc = c(0.5, 0.27), MethodCL = c("profile", 
                                                                                                                                                                       "Wald", "boot"), LevelCL = 0.95, roundstats = 3, cex = 0.75, 
                         cex.leg = 0.75, ncol.leg = 1, lty.overall = 1, lwd.overall = 2, 
                         col.overall = "red", RMprompt = F, xlab = "Estimated risk", 
                         ylab = "Observed proportion", xlim = c(0, 1), ylim = c(0, 
                                                                                1), d0lab = "0", d1lab = "1", cex.d01 = 0.7, dist.label = 0.04, 
                         line.bins = -0.05, dist.label2 = 0.03, las = 1, length.seg = 1, 
                         y.intersp = 1, lty.ideal = 1, col.ideal = "black", lwd.ideal = 1.75, 
                         lty.centers = NULL, lwd.centers = NULL, col.centers = NULL, 
                         EpsGrad = 0.001, fNonConv = c("warning", "stop"), Controlglmer = glmerControl(optimizer = "bobyqa"), 
                         ...) 
{
  Argz = as.list(match.call())[-1]
  CalibrLines = match.arg(CalibrLines)
  MethCL = match.arg(MethodCL)
  p = eval(Argz$p, data)
  y = eval(Argz$y, data)
  center = eval(Argz$center, data)
  LP = logit(p)
  if (is.factor(center)) 
    center = as.character(center)
  if (length(unique(center)) == 1) 
    stop("There is only one center, hence a random effects model should not be used here.")
  if (length(unique(center)) < 5) 
    warning("There are less than 5 centers, consider using a different method.", 
            immediate. = T)
  if (LogCal & flexible) 
    stop("LogCal and flexible cannot both be true.")
  Df = data.frame(y = y, p = p, LP = LP, center = center)
  Df = Df[with(Df, order(center, p)), ]
  RmCenter = dlply(Df, .(center), function(x) if (sum(x$y == 
                                                      1) < 10 | sum(x$y == 0) < 10) 
    unique(x$center)
    else NULL)
  RmCenter = unname(unlist(RmCenter))
  if (!is.null(RmCenter)) {
    if (RMprompt) {
      Cmmnd = readline(paste0("The center(s) ", paste0(RmCenter, 
                                                       collapse = ", "), " have less than 10 (non-)events and these will be", 
                              " removed. Do you want to continue? (y/n)   "))
      if (Cmmnd == "n") 
        stop("Function was stopped by the user.")
    }
    Df = Df[!(Df$center %in% RmCenter), ]
  }
  IncludedCenters = unique(Df$center)
  centers = IncludedCenters
  nCenters = length(IncludedCenters)
  if (CalibrLines != "overall") {
    if (all(sapply(list(lty.centers, lwd.centers, col.centers), 
                   is.null))) {
      lty.centers = rep(1, nCenters)
      lwd.centers = rep(1, nCenters)
      col.centers = seq_len(nCenters)
    }
    else if (sum(!sapply(list(lty.centers, lwd.centers, col.centers), 
                         is.null)) != 3) {
      stop("If you specify one of the arguments, please also specify the others.")
    }
    else {
      if (any(sapply(list(lty.centers, lwd.centers, col.centers), 
                     function(x) length(x) < nCenters))) 
        stop("The vector length of lty.centers, lwd.centers and col.centers is less than the number of centers.")
      FixL = function(x) x[1:nCenters]
      lty.centers = FixL(lty.centers)
      col.centers = FixL(col.centers)
      lwd.centers = FixL(lwd.centers)
    }
  }
  Df$center = factor(Df$center)
  contrasts(Df$center) = "contr.sum"
  LogMM = glmer(y ~ LP + (LP | center), data = Df, family = "binomial", 
                control = Controlglmer)
  CalSlope = fixef(LogMM)[2]
  cat("\n\nComputing confidence interval for calibration slope...\n\n")
  CalSlope = c(CalSlope, confint.merMod(LogMM, parm = "LP", 
                                        quiet = T, method = MethodCL, level = LevelCL))
  LogMM2 = glmer(y ~ 1 + (1 | center), data = Df, family = "binomial", 
                 offset = LP, control = Controlglmer)
  CalInterc = fixef(LogMM2)
  cat("\n\nComputing confidence interval for calibration intercept...\n\n")
  CalInterc = c(CalInterc, confint(LogMM2, parm = "(Intercept)", 
                                   quiet = T, method = MethodCL, level = LevelCL))

  
  if(CalibrLines != "overall"){
    par(mar = c(5,5,1,13), xpd=TRUE, pty = 's') 
  } else{
    par(pty = 's')
  }
  
  plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
       ylab = ylab, las = las, ...)
  clip(0, 1, 0, 1)
  abline(0, 1, lty = lty.ideal, col = col.ideal, lwd = lwd.ideal)
  do.call("clip", as.list(par()$usr))
  lt = lty.ideal
  lw.d = lwd.ideal
  all.col = col.ideal
  leg = "Ideal"
  NewX = data.frame(LP = seq(min(Df$LP), max(Df$LP), length = 500))
  p = inv.logit(NewX$LP)
  X = cbind(1, logit(p))
  FE = fixef(LogMM)
  if (CalibrLines != "overall") {
    PerC = data.frame(LP = rep(NewX$LP, length(IncludedCenters)), 
                      center = sort(rep(IncludedCenters, 500)))
    EstC = predict(LogMM, newdata = PerC, re.form = ~(LP | 
                                                        center), allow.new.levels = T, type = "response")
    ResultsC = cbind.data.frame(x = inv.logit(PerC$LP), y = EstC, 
                                center = PerC$center)
    ResultsC = split(ResultsC, ResultsC$center)
  }
  if (CalibrLines == "centers") {
    for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                              i)], lty = lty.centers[which(IncludedCenters == i)], 
                                     lwd = lwd.centers[which(IncludedCenters == i)])
    lt = c(lt, lty.centers)
    lw.d = c(lw.d, lwd.centers)
    all.col = c(all.col, col.centers)
    leg = c(leg, as.character(IncludedCenters))
  }
  else {
    p = inv.logit(NewX$LP)
    X = cbind(1, logit(p))
    FE = fixef(LogMM)
    if (!flexible) {
      OverallCal = inv.logit(X[order(p), ] %*% FE)
      p = p[order(p)]
      lines(p, OverallCal, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
    }
    else {
      Lfit = loess(y ~ p, Df)
      x = Lfit$x
      y = Lfit$fitted
      y = y[order(x)]
      x = x[order(x)]
      lines(x, y, lwd = lwd.overall, lty = lty.overall, 
            col = col.overall)
    }
    lt = c(lt, lty.overall)
    lw.d = c(lw.d, lwd.overall)
    all.col = c(all.col, col.overall)
    leg = c(leg, "Overall")
    if (CalibrLines == "both") {
      for (i in IncludedCenters) lines(ResultsC[[i]], col = col.centers[which(IncludedCenters == 
                                                                                i)], lty = lty.centers[which(IncludedCenters == 
                                                                                                               i)], lwd = lwd.centers[which(IncludedCenters == 
                                                                                                                                              i)])
      lt = c(lt, lty.centers)
      lw.d = c(lw.d, lwd.centers)
      all.col = c(all.col, col.centers)
      leg = c(leg, as.character(IncludedCenters))
    }
  }
  lp = legendloc
  lp = list(x = lp[1], y = lp[2])
  legend("topright", leg, lty = lt, cex = cex.leg, bty = "n", lwd = lw.d, 
         col = all.col, y.intersp = y.intersp, ncol = ncol.leg, inset = c(-.68, 0), xpd=NA) # Oorspronkelijk: c(-.39, 0)
  if (dostats) {
    if (CalibrLines == "centers"){
      stats.2 <- paste("")
      text(statloc[1], statloc[2], stats.2, pos = 4, cex = 0.75)
    } else{
      stats.2 <- matrix(ncol = 2, nrow = 2)
      colnames(stats.2) <- c("", "Estimate (95% CI)")
      stats.2[1, ] <- c("Intercept", paste0(format(round(CalInterc[1], 2), nsmall = 2), " (", format(round(CalInterc[2], 2), nsmall = 2), " to ", format(round(CalInterc[3], 2), nsmall = 2), ")"))
      stats.2[2, ] <- c("Slope", paste0(format(round(CalSlope[1], 2), nsmall = 2), " (", format(round(CalSlope[2], 2), nsmall = 2), " to ", format(round(CalSlope[3], 2), nsmall = 2), ")"))
      
      addtable2plot(x = statloc[1], y = statloc[2], table = stats.2, display.colnames = TRUE, cex = 0.75)
      
    }
  }
  
  Performance = rbind(CalInterc, CalSlope)
  rownames(Performance) = c("Calibration intercept", "Calibration slope"
                            )
  colnames(Performance) = c("Point estimate", "LCL", "UCL")
  return(structure(list(Performance = Performance,  
                        included = unique(Df$center), FitCalIntercept = LogMM2, 
                        FitCalSlope = LogMM, ConfLevel = LevelCL), class = "RE_ValProb"))
}

NB.imp <- function(pred, outcome,  data, imp,
                   sequence = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, imp = imp, stringsAsFactors = F)
  NRimp <- length(unique(Df$imp))
  
  
  ConfusionList <- list()
  for(i in 1:length(sequence)){
    ConfusionList[[i]] <- list()
  }
  
  ConfusionImp <- list()
  for(i in 1:NRimp){
    ConfusionImp[[i]] <- list()
  }
  
  for(k in 1:NRimp){
    #cat("Imputation", k, "of", NRimp)
    for(i in 1:length(sequence)){
      threshold <- sequence[i]
      
      Confusion <- matrix(nrow = NRimp, ncol = 9)
      Confusion <- data.frame(Confusion)
      colnames(Confusion) <- c('CutOff', 'TN', 'TP', 'FP', 'FN', 'cases', 'controls', 'NBTA', 'prevalence')
      Confusion$CutOff <- threshold
      
      
      CM <- confusion.matrix(obs = Df$y[Df$imp == k], pred = Df$p[Df$imp == k], threshold = threshold)
      Confusion$TN <- CM[1,1]
      Confusion$TP <- CM[2,2]
      Confusion$FP <- CM[2,1]
      Confusion$FN <- CM[1,2]
      
      Confusion$cases <- Confusion$TP + Confusion$FN
      Confusion$controls <- Confusion$TN + Confusion$FP
      
      Confusion$n <- Confusion$cases + Confusion$controls
      Confusion$NB <- Confusion$TP / Confusion$n - Confusion$FP / Confusion$n * (threshold / (1 - threshold))
      Confusion$prevalence <- Confusion$cases/Confusion$n
      Confusion$NBTA <- Confusion$prevalence - ((1-Confusion$prevalence) * (threshold/(1-threshold)))
      
      
      ConfusionList[[i]] <- Confusion
    }
    #cat("rbindlist toepassen")
    ConfusionImp[[k]] <- rbindlist(ConfusionList, fill = TRUE)
  }
  
  ConfusionLong <- rbindlist(ConfusionImp, fill = TRUE)
  ConfusionSum <- summaryBy(cbind(TP, TN, cases, controls, FP, FN, NB, NBTA) ~ cbind(CutOff), data = ConfusionLong, FUN = mean)
  
  return(structure(list(Results = ConfusionSum)))
}

PDIimp.centre.BAYES <- function(pred, outcome, center, imp, data, method.MA = "BAYES", titleGraph = "PDI per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  PDIimp <- list()
  for(i in 1:NRimp){
    PDIimp[[i]] <- list()
  }
  
  # PDI per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    PDIcenter <- matrix(ncol = 5, nrow = length(centers))
    colnames(PDIcenter) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
    PDIcenter <- data.frame(PDIcenter)
    PDIcenter$Center <- centers
    
    for(i in seq_along(centers)){
      PDIcenter[i,3 ] <- ests(y=Df$y[Df$center == centers[i] & Df$imp==j], d=Df[Df$center == centers[i] &Df$imp == j,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)[1]
      PDIcenter [i,4:5] <- ests(y=Df$y[Df$center == centers[i] & Df$imp==j], d=Df[Df$center == centers[i] &Df$imp == j,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)[3]
      PDIcenter[i, 2]   <- nrow(Df[Df$center == centers[i] & Df$imp == j,])
      
      
      PDIcenter$logit.PDI[i] <- logit(PDIcenter$PDI[i])
      PDIcenter$logit.se[i]  <- (logit(PDIcenter$PDI[i]) - logit(PDIcenter$LL[i]))/1.96
      PDIcenter$logit.var[i] <- PDIcenter$logit.se[i]^2
      
      
    }
    PDIcenter <- PDIcenter[order(PDIcenter$SampleSize, decreasing = TRUE),]
    
    PDIimp[[j]] <- PDIcenter
  }
  
  PDIimpLong <- rbindlist(PDIimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  PDIcombined <- matrix(ncol = 5, nrow = length(centers))
  colnames(PDIcombined) <- c('Center', 'SampleSize',  'logit.PDI', 'logit.LL', 'logit.UL')
  PDIcombined <- data.frame(PDIcombined)
  PDIcombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    PDIcombined$SampleSize[i] <- unique(PDIimpLong$SampleSize[PDIimpLong$Center == centers[i]])
    PDIcombined[i, 3] <- mean(PDIimpLong$logit.PDI[PDIimpLong$Center == centers[i]])
    WithinVar <- mean(PDIimpLong$logit.var[PDIimpLong$Center == centers[i]])
    BetweenVar <- var(PDIimpLong$logit.PDI[PDIimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    PDIcombined$PooledSE[i] <- sqrt(PooledVar)
    PDIcombined$logit.LL[i] <- PDIcombined$logit.PDI[i] - 1.96*PDIcombined$PooledSE[i]
    PDIcombined$logit.UL[i] <- PDIcombined$logit.PDI[i] + 1.96*PDIcombined$PooledSE[i]
  }
  
  PDIcombined$PDI <- inv.logit(PDIcombined$logit.PDI)
  PDIcombined$LL <- inv.logit(PDIcombined$logit.LL)
  PDIcombined$UL <- inv.logit(PDIcombined$logit.UL)
  PDIcombined <- PDIcombined[order(PDIcombined$SampleSize, decreasing = TRUE),]
  
  PDIoverall <- matrix(nrow = 2, ncol = 5)
  colnames(PDIoverall) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  PDIoverall <- data.frame(PDIoverall)
  PDIoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  PDIoverall$SampleSize <- nrow(Df[Df$imp == 1,])
  
  # Meta-analyse voor overall estimate
  fit.RE = uvmeta(PDIcombined$logit.PDI, r.se = PDIcombined$PooledSE, method = method.MA)
  
  
  PDIoverall$PDI[1] <- inv.logit(as.numeric(fit.RE[7]))
  PDIoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  PDIoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  PDIoverall$PDI[2] <- inv.logit(as.numeric(fit.RE[7]))
  PDIoverall$LL[2] <- inv.logit(as.numeric(fit.RE[13]))
  PDIoverall$UL[2] <- inv.logit(as.numeric(fit.RE[14]))
  PDIoverall
  
  NAforest <- matrix(nrow = 1, ncol = 5)
  colnames(NAforest) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  PDI <- rbind(PDIcombined[, c('Center', 'SampleSize',  'PDI', 'LL', 'UL')], NAforest, NAforest, PDIoverall)
  
  # Layout for forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(PDI)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- PDI$Center[i]
  }
  RRpdi <- 1:nrobs
  for(i in 1:nrobs){
    RRpdi[i] <- paste(format(round(PDI$PDI[i], 2), nsmall = 2), " (", format(round(PDI$LL[i], 2), nsmall = 2), " to ", format(round(PDI$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- PDI$SampleSize[i]
  }
  
  Labels <- c('Centre', 'PDI (95% CI)', 'N') #(prev)')
  Combined <- cbind(RRcenter, RRpdi, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  PDIna <- rbind(NAforest, NAforest, PDI)
  
  return(structure(list(Performance = PDIoverall, ModelFit = fit.RE, PDIcenters = PDIcombined[, c('Center', 'SampleSize',  'PDI', 'LL', 'UL')], IncludedCenters = centers, dataPlot = PDIna, Plot = Tabletext, data = Df)))
  
}



PDI.centre.BAYES <- function(pred, outcome, center,data, method.MA = "BAYES", titleGraph = "PDI per center"){
  
  arguments <- as.list(match.call())[-1]
  #pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  
  
  Df = data.frame(p = pred, y = outcome, center = center,  stringsAsFactors = F)
  centers <- unique(Df$center)
  
  
  
  # PDI per center
  
  
  PDIcenter <- matrix(ncol = 5, nrow = length(centers))
  colnames(PDIcenter) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  PDIcenter <- data.frame(PDIcenter)
  PDIcenter$Center <- centers
  
  for(i in seq_along(centers)){
    PDIcenter[i,3 ] <- ests(y=Df$y[Df$center == centers[i] ], d=Df[Df$center == centers[i] ,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)[1]
    PDIcenter [i,4] <- ests(y=Df$y[Df$center == centers[i] ], d=Df[Df$center == centers[i] ,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)$interval[1]
    PDIcenter [i,5] <- ests(y=Df$y[Df$center == centers[i] ], d=Df[Df$center == centers[i] ,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)$interval[2]
    PDIcenter[i, 2]   <- nrow(Df[Df$center == centers[i] ,])
    
    
    PDIcenter$logit.PDI[i] <- logit(PDIcenter$PDI[i])
    PDIcenter$logit.se[i]  <- (logit(PDIcenter$PDI[i]) - logit(PDIcenter$LL[i]))/1.96
    PDIcenter$logit.var[i] <- PDIcenter$logit.se[i]^2
    
    
  }
  
  PDIcenter$logit.PDI <- logit(PDIcenter$PDI)
  PDIcenter$logit.se <- (logit(PDIcenter$PDI)-logit(PDIcenter$LL))/1.96
  PDIcenter <- PDIcenter[order(PDIcenter$SampleSize, decreasing = TRUE),]
  
  
  PDIoverall <- matrix(nrow = 2, ncol = 5)
  colnames(PDIoverall) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  PDIoverall <- data.frame(PDIoverall)
  PDIoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  PDIoverall$SampleSize <- nrow(Df)
  
  
  # Meta-analyse voor overall estimate
  fit.RE = uvmeta(PDIcenter$logit.PDI, r.se = PDIcenter$logit.se, method = method.MA)
  
  
  PDIoverall$PDI[1] <- inv.logit(as.numeric(fit.RE[7]))
  PDIoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  PDIoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  PDIoverall$PDI[2] <- inv.logit(as.numeric(fit.RE[7]))
  PDIoverall$LL[2] <- inv.logit(as.numeric(fit.RE[13]))
  PDIoverall$UL[2] <- inv.logit(as.numeric(fit.RE[14]))
  PDIoverall
  
  NAforest <- matrix(nrow = 1, ncol = 5)
  colnames(NAforest) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  PDI <- rbind(PDIcenter[, c('Center', 'SampleSize',  'PDI', 'LL', 'UL')], NAforest,  PDIoverall)
  
  # Layout for forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(PDI)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- PDI$Center[i]
  }
  RRpdi <- 1:nrobs
  for(i in 1:nrobs){
    RRpdi[i] <- paste(format(round(PDI$PDI[i], 2), nsmall = 2), " (", format(round(PDI$LL[i], 2), nsmall = 2), " to ", format(round(PDI$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- PDI$SampleSize[i]
  }
  
  Labels <- c('Centre', 'PDI (95% CI)', 'N') #(prev)')
  Combined <- cbind(RRcenter, RRpdi, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  PDIna <- rbind(NAforest, NAforest, PDI)
  
  return(structure(list(Performance = PDIoverall, ModelFit = fit.RE, PDIcenters = PDIcenter, IncludedCenters = centers, dataPlot = PDIna, Plot = Tabletext, data = Df)))
  
}

PDI.centre <- function(pred, outcome, center,data, method.MA = "REML", titleGraph = "PDI per center"){
  
  arguments <- as.list(match.call())[-1]
  #pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  
  
  Df = data.frame(p = pred, y = outcome, center = center,  stringsAsFactors = F)
  centers <- unique(Df$center)
  
  
  
  # PDI per center
  
  
  PDIcenter <- matrix(ncol = 5, nrow = length(centers))
  colnames(PDIcenter) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  PDIcenter <- data.frame(PDIcenter)
  PDIcenter$Center <- centers
  
  for(i in seq_along(centers)){
    PDIcenter[i,3 ] <- ests(y=Df$y[Df$center == centers[i] ], d=Df[Df$center == centers[i] ,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)[1]
    PDIcenter [i,4] <- ests(y=Df$y[Df$center == centers[i] ], d=Df[Df$center == centers[i] ,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)$interval[1]
    PDIcenter [i,5] <- ests(y=Df$y[Df$center == centers[i] ], d=Df[Df$center == centers[i] ,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)$interval[2]
    PDIcenter[i, 2]   <- nrow(Df[Df$center == centers[i] ,])
    
    
    PDIcenter$logit.PDI[i] <- logit(PDIcenter$PDI[i])
    PDIcenter$logit.se[i]  <- (logit(PDIcenter$PDI[i]) - logit(PDIcenter$LL[i]))/1.96
    PDIcenter$logit.var[i] <- PDIcenter$logit.se[i]^2
    
    
  }
  
  PDIcenter$logit.PDI <- logit(PDIcenter$PDI)
  PDIcenter$logit.se <- (logit(PDIcenter$PDI)-logit(PDIcenter$LL))/1.96
  PDIcenter <- PDIcenter[order(PDIcenter$SampleSize, decreasing = TRUE),]
  
  
  PDIoverall <- matrix(nrow = 2, ncol = 5)
  colnames(PDIoverall) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  PDIoverall <- data.frame(PDIoverall)
  PDIoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  PDIoverall$SampleSize <- nrow(Df)
  
  
  # Meta-analyse voor overall estimate

  
  fit.RE = rma.uni(PDIcenter$logit.PDI, sei = PDIcenter$logit.se, method = method.MA)
  PI = predict(fit.RE, transf = transf.ilogit)
  
  PDIoverall$PDI[1] <- inv.logit(coef(fit.RE))
  PDIoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  PDIoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  PDIoverall$PDI[2] <- PI$pred
  PDIoverall$LL[2] <- PI$cr.lb
  PDIoverall$UL[2] <- PI$cr.ub
  PDIoverall
  
  NAforest <- matrix(nrow = 1, ncol = 5)
  colnames(NAforest) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  PDI <- rbind(PDIcenter[, c('Center', 'SampleSize',  'PDI', 'LL', 'UL')], NAforest,  PDIoverall)
  
  # Layout for forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(PDI)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- PDI$Center[i]
  }
  RRpdi <- 1:nrobs
  for(i in 1:nrobs){
    RRpdi[i] <- paste(format(round(PDI$PDI[i], 2), nsmall = 2), " (", format(round(PDI$LL[i], 2), nsmall = 2), " to ", format(round(PDI$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- PDI$SampleSize[i]
  }
  
  Labels <- c('Centre', 'PDI (95% CI)', 'N') #(prev)')
  Combined <- cbind(RRcenter, RRpdi, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  PDIna <- rbind(NAforest, NAforest, PDI)
  
  return(structure(list(Performance = PDIoverall, ModelFit = fit.RE, PDIcenters = PDIcenter, IncludedCenters = centers, dataPlot = PDIna, Plot = Tabletext, data = Df)))
  
}


AUC.IOTA <- function(pred, outcome, center, data, method.MA = "REML", titleGraph = "AUC per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  AUCcenter <- matrix(ncol = 6, nrow = length(centers))
  colnames(AUCcenter) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCcenter <- data.frame(AUCcenter)
  AUCcenter$Center <- centers
  
  # AUC per center
  for(i in seq_along(centers)){
    AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 1], Df$p[Df$center == centers[i] & Df$y == 3], method = "pepe")
    AUCcenter[i, 2]   <- nrow(Df[Df$center == centers[i],])
    AUCcenter[i, 3]   <- round(nrow(Df[Df$y == 3 & Df$center == centers[i],])/nrow(Df[Df$center == centers[i],])*100)
    
    
      AUCcenter$logit.AUC[i] <- logit(AUCcenter$AUC[i])
      AUCcenter$logit.se[i]  <- (logit(AUCcenter$AUC[i]) - logit(AUCcenter$LL[i]))/1.96
      AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2

  }
  
  
  AUCcenter$logit.AUC <- logit(AUCcenter$AUC)
  AUCcenter$logit.se  <- (logit(AUCcenter$AUC) - logit(AUCcenter$LL))/1.96
  AUCcenter <- AUCcenter[order(AUCcenter$SampleSize, decreasing = TRUE),]
  
  AUCoverall <- matrix(nrow = 2, ncol = 6)
  colnames(AUCoverall) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCoverall <- data.frame(AUCoverall)
  AUCoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  AUCoverall$SampleSize <- nrow(Df)
  AUCoverall$Prevalence <- round(nrow(Df[Df$y == 2,])/nrow(Df)*100)
  
  # Meta-analyse voor overall estimate
  fit.RE = rma.uni(AUCcenter$logit.AUC, sei = AUCcenter$logit.se, method = method.MA)
  PI = predict(fit.RE, transf = transf.ilogit)
  
  AUCoverall$AUC[1] <- inv.logit(coef(fit.RE))
  AUCoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  AUCoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  AUCoverall$AUC[2] <- PI$pred
  AUCoverall$LL[2] <- PI$cr.lb
  AUCoverall$UL[2] <- PI$cr.ub
  AUCoverall
  
  NAforest <- matrix(nrow = 1, ncol = 6)
  colnames(NAforest) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  AUC <- rbind(AUCcenter[, 1:6], NAforest, AUCoverall)
  
  # Layout forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(AUC)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- AUC$Center[i]
  }
  RRauc <- 1:nrobs
  for(i in 1:nrobs){
    RRauc[i] <- paste(format(round(AUC$AUC[i], 2), nsmall = 2), " (", format(round(AUC$LL[i], 2), nsmall = 2), " to ", format(round(AUC$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- AUC$SampleSize[i]
  }
  
  Labels <- c('Centre', 'AUC (95% CI)', 'N')
  Combined <- cbind(RRcenter, RRauc, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  AUCna <- rbind(NAforest, NAforest, AUC)
  
  return(structure(list(Performance = AUCoverall, ModelFit = fit.RE, AUCcenters = AUCcenter, IncludedCenters = centers, dataPlot = AUCna, Plot = Tabletext, data = Df)))
  
}

AUC.IOTA.BAYES <- function(pred, outcome, center, data, method.MA = "BAYES", titleGraph = "AUC per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  AUCcenter <- matrix(ncol = 6, nrow = length(centers))
  colnames(AUCcenter) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCcenter <- data.frame(AUCcenter)
  AUCcenter$Center <- centers
  
  # AUC per center
  for(i in seq_along(centers)){
    AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 4], Df$p[Df$center == centers[i] & Df$y == 5], method = "pepe")
    AUCcenter[i, 2]   <- nrow(Df[Df$center == centers[i],])
    AUCcenter[i, 3]   <- round(nrow(Df[Df$y == 4 & Df$center == centers[i],])/nrow(Df[Df$center == centers[i],])*100)
    
    
    AUCcenter$logit.AUC[i] <- logit(AUCcenter$AUC[i])
    AUCcenter$logit.se[i]  <- (logit(AUCcenter$AUC[i]) - logit(AUCcenter$LL[i]))/1.96
    AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
    
  }
  
  
  AUCcenter$logit.AUC <- logit(AUCcenter$AUC)
  AUCcenter$logit.se  <- (logit(AUCcenter$AUC) - logit(AUCcenter$LL))/1.96
  AUCcenter <- AUCcenter[order(AUCcenter$SampleSize, decreasing = TRUE),]
  
  AUCoverall <- matrix(nrow = 2, ncol = 6)
  colnames(AUCoverall) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCoverall <- data.frame(AUCoverall)
  AUCoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  AUCoverall$SampleSize <- nrow(Df)
  AUCoverall$Prevalence <- round(nrow(Df[Df$y == 2,])/nrow(Df)*100)
  
  # Meta-analyse voor overall estimate
  # Meta-analyse voor overall estimate
  fit.RE = uvmeta(AUCcenter$logit.AUC, r.se = AUCcenter$logit.se, method = method.MA)
  
  
  AUCoverall$AUC[1] <- inv.logit(as.numeric(fit.RE[7]))
  AUCoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  AUCoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  AUCoverall$AUC[2] <- inv.logit(as.numeric(fit.RE[7]))
  AUCoverall$LL[2] <- inv.logit(as.numeric(fit.RE[13]))
  AUCoverall$UL[2] <- inv.logit(as.numeric(fit.RE[14]))
  AUCoverall
  

  
  NAforest <- matrix(nrow = 1, ncol = 6)
  colnames(NAforest) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  AUC <- rbind(AUCcenter[, 1:6], NAforest, AUCoverall)
  
  # Layout forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(AUC)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- AUC$Center[i]
  }
  RRauc <- 1:nrobs
  for(i in 1:nrobs){
    RRauc[i] <- paste(format(round(AUC$AUC[i], 2), nsmall = 2), " (", format(round(AUC$LL[i], 2), nsmall = 2), " to ", format(round(AUC$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- AUC$SampleSize[i]
  }
  
  Labels <- c('Centre', 'AUROC (95% CI)', 'N')
  Combined <- cbind(RRcenter, RRauc, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  AUCna <- rbind(NAforest, NAforest, AUC)
  
  return(structure(list(Performance = AUCoverall, ModelFit = fit.RE, AUCcenters = AUCcenter, IncludedCenters = centers, dataPlot = AUCna, Plot = Tabletext, data = Df)))
  
}

PDIimp.centre <- function(pred, outcome, center, imp, data, method.MA = "REML", titleGraph = "PDI per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  PDIimp <- list()
  for(i in 1:NRimp){
    PDIimp[[i]] <- list()
  }
  
  # PDI per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    PDIcenter <- matrix(ncol = 5, nrow = length(centers))
    colnames(PDIcenter) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
    PDIcenter <- data.frame(PDIcenter)
    PDIcenter$Center <- centers
    
    for(i in seq_along(centers)){
      PDIcenter[i,3 ] <- ests(y=Df$y[Df$center == centers[i] & Df$imp==j], d=Df[Df$center == centers[i] &Df$imp == j,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)[1]
      PDIcenter [i,4:5] <- ests(y=Df$y[Df$center == centers[i] & Df$imp==j], d=Df[Df$center == centers[i] &Df$imp == j,1:5], acc="pdi", level=0.95, method="prob", k=5, B=5)[3]
      PDIcenter[i, 2]   <- nrow(Df[Df$center == centers[i] & Df$imp == j,])
      
      
      PDIcenter$logit.PDI[i] <- logit(PDIcenter$PDI[i])
      PDIcenter$logit.se[i]  <- (logit(PDIcenter$PDI[i]) - logit(PDIcenter$LL[i]))/1.96
      PDIcenter$logit.var[i] <- PDIcenter$logit.se[i]^2
      
      
    }
    PDIcenter <- PDIcenter[order(PDIcenter$SampleSize, decreasing = TRUE),]
    
    PDIimp[[j]] <- PDIcenter
  }
  
  PDIimpLong <- rbindlist(PDIimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  PDIcombined <- matrix(ncol = 5, nrow = length(centers))
  colnames(PDIcombined) <- c('Center', 'SampleSize',  'logit.PDI', 'logit.LL', 'logit.UL')
  PDIcombined <- data.frame(PDIcombined)
  PDIcombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    PDIcombined$SampleSize[i] <- unique(PDIimpLong$SampleSize[PDIimpLong$Center == centers[i]])
    PDIcombined[i, 3] <- mean(PDIimpLong$logit.PDI[PDIimpLong$Center == centers[i]])
    WithinVar <- mean(PDIimpLong$logit.var[PDIimpLong$Center == centers[i]])
    BetweenVar <- var(PDIimpLong$logit.PDI[PDIimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    PDIcombined$PooledSE[i] <- sqrt(PooledVar)
    PDIcombined$logit.LL[i] <- PDIcombined$logit.PDI[i] - 1.96*PDIcombined$PooledSE[i]
    PDIcombined$logit.UL[i] <- PDIcombined$logit.PDI[i] + 1.96*PDIcombined$PooledSE[i]
  }
  
  PDIcombined$PDI <- inv.logit(PDIcombined$logit.PDI)
  PDIcombined$LL <- inv.logit(PDIcombined$logit.LL)
  PDIcombined$UL <- inv.logit(PDIcombined$logit.UL)
  PDIcombined <- PDIcombined[order(PDIcombined$SampleSize, decreasing = TRUE),]
  
  PDIoverall <- matrix(nrow = 2, ncol = 5)
  colnames(PDIoverall) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  PDIoverall <- data.frame(PDIoverall)
  PDIoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  PDIoverall$SampleSize <- nrow(Df[Df$imp == 1,])
  
  # Meta-analyse voor overall estimate
  fit.RE = rma.uni(PDIcombined$logit.PDI, sei = PDIcombined$PooledSE, method = "REML")
  PI = predict(fit.RE, transf = transf.ilogit)
  
  PDIoverall$PDI[1] <- inv.logit(coef(fit.RE))
  PDIoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  PDIoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  PDIoverall$PDI[2] <- PI$pred
  PDIoverall$LL[2] <- PI$cr.lb
  PDIoverall$UL[2] <- PI$cr.ub
  PDIoverall
  
  NAforest <- matrix(nrow = 1, ncol = 5)
  colnames(NAforest) <- c('Center', 'SampleSize',  'PDI', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  PDI <- rbind(PDIcombined[, c('Center', 'SampleSize',  'PDI', 'LL', 'UL')], NAforest, NAforest, PDIoverall)
  
  # Layout for forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(PDI)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- PDI$Center[i]
  }
  RRpdi <- 1:nrobs
  for(i in 1:nrobs){
    RRpdi[i] <- paste(format(round(PDI$PDI[i], 2), nsmall = 2), " (", format(round(PDI$LL[i], 2), nsmall = 2), " to ", format(round(PDI$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- PDI$SampleSize[i]
  }
  
  Labels <- c('Centre', 'PDI (95% CI)', 'N') #(prev)')
  Combined <- cbind(RRcenter, RRpdi, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  PDIna <- rbind(NAforest, NAforest, PDI)
  
  return(structure(list(Performance = PDIoverall, ModelFit = fit.RE, PDIcenters = PDIcombined[, c('Center', 'SampleSize',  'PDI', 'LL', 'UL')], IncludedCenters = centers, dataPlot = PDIna, Plot = Tabletext, data = Df)))
  
}


AUCimp.IOTA.BAYES <- function(pred, outcome, center, imp, data, method.MA = "BAYES", titleGraph = "AUC per center"){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  NRimp <- length(unique(Df$imp))
  
  AUCimp <- list()
  for(i in 1:NRimp){
    AUCimp[[i]] <- list()
  }
  
  PrevalenceOverall <- matrix(nrow = NRimp, ncol = 1)
  
  # AUC per center
  for(j in 1:NRimp){
    cat("Imputation", j, " of ", NRimp, "\n\n")
    
    AUCcenter <- matrix(ncol = 6, nrow = length(centers))
    colnames(AUCcenter) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
    AUCcenter <- data.frame(AUCcenter)
    AUCcenter$Center <- centers
    
    for(i in seq_along(centers)){
      AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 4 & Df$imp == j], Df$p[Df$center == centers[i] & Df$y ==5 & Df$imp == j], method = "pepe")
      AUCcenter[i, 2]   <- nrow(Df[Df$center == centers[i] & Df$imp == j,])
      AUCcenter[i, 3]   <- round(nrow(Df[Df$y == 4 & Df$center == centers[i] & Df$imp == j,])/nrow(Df[Df$center == centers[i] & Df$imp == j,])*100)
      
      ## Additional part for AUCs of 1
      if(AUCcenter[i, 4] == 1){
        AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 4 & Df$imp == j], Df$p[Df$center == centers[i] & Df$y == 5 & Df$imp == j], method = "newcombe") # Newcombe ipv pepe
      } else{
        AUCcenter[i, 4:6] <- auc.nonpara.mw(Df$p[Df$center == centers[i] & Df$y == 4 & Df$imp == j], Df$p[Df$center == centers[i] & Df$y == 5 & Df$imp == j], method = "pepe")
      }
      
      if(AUCcenter$AUC[i] != 1){
        AUCcenter$logit.AUC[i] <- logit(AUCcenter$AUC[i])
        AUCcenter$logit.se[i]  <- (logit(AUCcenter$AUC[i]) - logit(AUCcenter$LL[i]))/1.96
        AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
      } else{
        AUCcenter$logit.AUC[i] <- logit(0.999)
        AUCcenter$logit.se[i]  <- (logit(0.999) - logit(AUCcenter$LL[i]))/1.96
        AUCcenter$logit.var[i] <- AUCcenter$logit.se[i]^2
      }
      
    }
    AUCcenter <- AUCcenter[order(AUCcenter$SampleSize, decreasing = TRUE),]
    
    AUCimp[[j]] <- AUCcenter
    
    PrevalenceOverall[j] <- round(nrow(Df[Df$y == 4 & Df$imp == j,])/nrow(Df[Df$imp == j,])*100)
  }
  
  AUCimpLong <- rbindlist(AUCimp, fill = TRUE)
  
  # Combine results with Rubin's rule
  AUCcombined <- matrix(ncol = 6, nrow = length(centers))
  colnames(AUCcombined) <- c('Center', 'SampleSize', 'Prevalence', 'logit.AUC', 'logit.LL', 'logit.UL')
  AUCcombined <- data.frame(AUCcombined)
  AUCcombined$Center <- centers
  
  
  for(i in seq_along(centers)){
    AUCcombined$SampleSize[i] <- unique(AUCimpLong$SampleSize[AUCimpLong$Center == centers[i]])
    AUCcombined$Prevalence[i] <- round(mean(AUCimpLong$Prevalence[AUCimpLong$Center == centers[i]]))
    AUCcombined[i, 4] <- mean(AUCimpLong$logit.AUC[AUCimpLong$Center == centers[i]])
    WithinVar <- mean(AUCimpLong$logit.var[AUCimpLong$Center == centers[i]])
    BetweenVar <- var(AUCimpLong$logit.AUC[AUCimpLong$Center == centers[i]])
    PooledVar <- WithinVar + BetweenVar + BetweenVar/NRimp
    AUCcombined$PooledSE[i] <- sqrt(PooledVar)
    AUCcombined$logit.LL[i] <- AUCcombined$logit.AUC[i] - 1.96*AUCcombined$PooledSE[i]
    AUCcombined$logit.UL[i] <- AUCcombined$logit.AUC[i] + 1.96*AUCcombined$PooledSE[i]
  }
  
  AUCcombined$AUC <- inv.logit(AUCcombined$logit.AUC)
  AUCcombined$LL <- inv.logit(AUCcombined$logit.LL)
  AUCcombined$UL <- inv.logit(AUCcombined$logit.UL)
  AUCcombined <- AUCcombined[order(AUCcombined$SampleSize, decreasing = TRUE),]
  
  AUCoverall <- matrix(nrow = 2, ncol = 6)
  colnames(AUCoverall) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  AUCoverall <- data.frame(AUCoverall)
  AUCoverall$Center <- c("Meta-analysis", "95% Prediction interval")
  AUCoverall$SampleSize <- nrow(Df[Df$imp == 1,])
  AUCoverall$Prevalence <- round(mean(PrevalenceOverall))
  
  # Meta-analyse voor overall estimate
  
  fit.RE=uvmeta(r=AUCcombined$logit.AUC, r.se=AUCcombined$PooledSE, method=method.MA)
  
  
  AUCoverall$AUC[1] <- inv.logit(as.numeric(fit.RE[7]))
  AUCoverall$LL[1] <- inv.logit(fit.RE$ci.lb)
  AUCoverall$UL[1] <- inv.logit(fit.RE$ci.ub)
  AUCoverall$AUC[2] <- inv.logit(as.numeric(fit.RE[7]))
  AUCoverall$LL[2] <- inv.logit(as.numeric(fit.RE[13]))
  AUCoverall$UL[2] <- inv.logit(as.numeric(fit.RE[14]))
  AUCoverall
  
  
  NAforest <- matrix(nrow = 1, ncol = 6)
  colnames(NAforest) <- c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')
  NAforest <- data.frame(NAforest)
  NAforest <- NA
  
  AUC <- rbind(AUCcombined[, c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')], NAforest, NAforest, AUCoverall)
  
  # Layout for forestplot
  RRcenter <- c('Other', 'Meta-analysis')
  nrobs <- nrow(AUC)
  RRcenter <- 1:nrobs
  for(i in 1:nrobs){
    RRcenter[i] <- AUC$Center[i]
  }
  RRauc <- 1:nrobs
  for(i in 1:nrobs){
    RRauc[i] <- paste(format(round(AUC$AUC[i], 2), nsmall = 2), " (", format(round(AUC$LL[i], 2), nsmall = 2), " to ", format(round(AUC$UL[i], 2), nsmall = 2), ")", sep = "")
  }
  RRprev <- 1:nrobs
  for(i in 1:nrobs){
    RRprev[i] <- AUC$SampleSize[i]
  }
  
  Labels <- c('Centre', 'AUROC (95% CI)', 'N') #(prev)')
  Combined <- cbind(RRcenter, RRauc, RRprev)
  Tabletext <- rbind(Labels, NAforest, Combined)
  
  AUCna <- rbind(NAforest, NAforest, AUC)
  
  return(structure(list(Performance = AUCoverall, ModelFit = fit.RE, AUCcenters = AUCcombined[, c('Center', 'SampleSize', 'Prevalence', 'AUC', 'LL', 'UL')], IncludedCenters = centers, dataPlot = AUCna, Plot = Tabletext, data = Df)))
  
}



## Without multiple imputation
DataWinBugs <- function(pred, outcome, center, data, 
                        sequence = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, stringsAsFactors = F)
  
  centers <- unique(Df$center)
  
  ConfusionList <- list()
  
  for(i in 1:length(sequence)){
    ConfusionList[[i]] <- list()
  }
  
  
  for(i in 1:length(sequence)){
    threshold <- sequence[i]
    
    Confusion <- matrix(nrow = length(centers), ncol = 8)
    Confusion <- data.frame(Confusion)
    colnames(Confusion) <- c('Center', 'CutOff', 'TN', 'TP', 'FP', 'FN', 'cases', 'controls')
    Confusion$CutOff <- threshold
    
    for(j in seq_along(centers)){
      
      Confusion$Center[j] <- centers[j]
      
      predictedvalues <- ifelse(Df$p[Df$center == centers[j]]>threshold,1,0)
      CM <- table(predictedvalues, Df$y[Df$center == centers[j]])
      
      #CM <- confusion.matrix(obs = Df$y[Df$center == centers[j]], pred = Df$p[Df$center == centers[j]], threshold = threshold)
      Confusion$TN[j] <- CM[1,1]
      Confusion$TP[j] <- CM[2,2]
      Confusion$FP[j] <- CM[2,1]
      Confusion$FN[j] <- CM[1,2]
      
      Confusion$cases <- Confusion$TP + Confusion$FN
      Confusion$controls <- Confusion$TN + Confusion$FP
      
      Confusion$n <- Confusion$cases + Confusion$controls
      Confusion$NB <- Confusion$TP / Confusion$n - Confusion$FP / Confusion$n * (threshold / (1 - threshold))
      
    }
    
    ConfusionList[[i]] <- Confusion
    
    
  }
  return(structure(list(Results = ConfusionList)))
}

#clinical utility
DataWinBugs.imp.CM <- function(pred, outcome, center, data, imp,
                            sequence = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  center <- eval(arguments$center, data)
  imp <- eval(arguments$imp, data)
  
  Df = data.frame(p = pred, y = outcome, center = center, imp = imp, stringsAsFactors = F)
  NRimp <- length(unique(Df$imp))
  
  centers <- unique(Df$center)
  
  ConfusionList <- list()
  for(i in 1:length(sequence)){
    ConfusionList[[i]] <- list()
  }
  
  ConfusionImp <- list()
  for(i in 1:NRimp){
    ConfusionImp[[i]] <- list()
  }
  
  for(k in 1:NRimp){
    cat("Imputation", k, "of", NRimp)
    for(i in 1:length(sequence)){
      threshold <- sequence[i]
      
      Confusion <- matrix(nrow = length(centers), ncol = 8)
      Confusion <- data.frame(Confusion)
      colnames(Confusion) <- c('Center', 'CutOff', 'TN', 'TP', 'FP', 'FN', 'cases', 'controls')
      Confusion$CutOff <- threshold
      
      for(j in seq_along(centers)){
        
        Confusion$Center[j] <- centers[j]
        
        predictedvalues <- ifelse(Df$p[Df$center == centers[j] & Df$imp == k]>threshold,1,0)
        CM <- table(predictedvalues, Df$y[Df$center == centers[j] & Df$imp == k])
        
        #CM <- confusion.matrix(obs = Df$y[Df$center == centers[j] & Df$imp == k], pred = Df$p[Df$center == centers[j] & Df$imp == k], threshold = threshold)
        Confusion$TN[j] <- CM[1,1]
        Confusion$TP[j] <- CM[2,2]
        Confusion$FP[j] <- CM[2,1]
        Confusion$FN[j] <- CM[1,2]
        
        Confusion$cases <- Confusion$TP + Confusion$FN
        Confusion$controls <- Confusion$TN + Confusion$FP
        
        Confusion$n <- Confusion$cases + Confusion$controls
        Confusion$NB <- Confusion$TP / Confusion$n - Confusion$FP / Confusion$n * (threshold / (1 - threshold))
        
      }
      
      ConfusionList[[i]] <- Confusion
    }
    cat("rbindlist toepassen")
    ConfusionImp[[k]] <- rbindlist(ConfusionList, fill = TRUE)
  }
  
  ConfusionLong <- rbindlist(ConfusionImp, fill = TRUE)
  ConfusionSum <- summaryBy(cbind(TP, TN, cases, controls, FP, FN, NB) ~ cbind(CutOff, Center), data = ConfusionLong, FUN = mean)
  
  return(structure(list(Results = ConfusionSum)))
}
