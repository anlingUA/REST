source("RAD.R")
require(quadprog)

REST <- function(data_target, data_suspect,nbootstrap=1000){
  # Performs REST source selection for microbial source tracking for rarefied abundance counts.
  # Args:
  #   data_target:  array containing a single test/sink sample.
  #   data_suspect: data frame containing all control/source samples.
  #   nbootstrap:   the number of bootstrap resamples to generate for
  #                 calculating boostrap estimates and test statistics.
  # Returns:
  #   proportions of suspected profiles
  resultsW2 = RAD(data_target,data_suspect,n=nbootstrap)
  RAD_indicator = c(resultsW2[[2]][7,],0)
  data_RAD <- data_suspect[,RAD_indicator ==1] # pick RAD selected suspect data set to calculate med in next step
  med <- mean(sum(data_RAD))# rarefy library size  #median selected column sum
  # weighted linear regression y is evidence 
  data2 <- cbind(data_target,data_RAD) #combine data for regression
  unknown<-missing(data_target,data,med) 
  unknown <- as.data.frame(unknown)
  suspects<-cbind(data_RAD, unknown) #suspects with unknown 
  # stage 4
  proposed <- t(constrained(data_target,suspects))
  colnames(proposed) <- c(colnames(data_RAD),"unknown")  
  return(proposed)
}

missing <- function(data_target,data,med){
  Y <- data_target
  fit1 <- lm(data_target ~ ., data = data)
  e1 <- resid(fit1)
  fity <- fitted(fit1)
  fit2 <- lm(abs(e1)~fity)
  fitvalue <- fitted(fit2)
  wi <- 1/(fitvalue)^2
  fit3 <- lm(data_target ~ .,data = data, weights=wi)
  e3 <- resid(fit3)
  res <- e3 * sqrt(wi)
  min <- min(res)
  di <- res-min
  unknown <- (di*med)/sum(di)
  return(unknown)
}

constrained <- function(data_target,suspects){
  X <- as.matrix(suspects)
  Q <- solve(chol(t(X) %*% X))
  numberofsusp <- ncol(suspects) # how many suspects in total
  C <- cbind(rep(1,numberofsusp),diag(numberofsusp)) 
  b <- c(1,rep(0,numberofsusp))  
  EV <- data_target
  dvec = t(EV)%*%X
  constrained<-solve.QP( Q, dvec, C, b,meq = 1, factorized = TRUE)$solution
  return(constrained) 
}

### optional ICC

library(psych)
# data_in: put unknown (obtained from missing function above) from multiple replicates together; 
#          cols are replicates, rows are OTUs ;
# ICC(data_in, lmer=FALSE)
