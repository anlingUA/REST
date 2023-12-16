 # File: RAD.R
# Author: Kyle Carter
# Contact: kcarter@math.arizona.edu
# Contact: anling@arizona.edu

# This file contains utility functions for performing bootstrap resampling
# and generating RAD microbial source tracking for data with rarefied
# integer abundance counts.

#

require(robCompositions)

ZeroFilter=function(data){
  # Removes rows from an integer counts data frame that contain zero counts 
  # for all observations.
  #
  # Args:
  #   data:   a data frame containing integer count data with observations
  #           along columns and variables along rows.
  # Returns:
  #   A data frame containing variables with at least one positive count.
  
  name=names(data)
  filter=as.matrix(data)
  filter=filter[which(rowSums(filter)>0),]
  filter=data.frame(filter)
  names(filter)=name
  return(filter)
}


RMultiMix <- function(n, size, prob, w = NULL) {
  # Generates mixed multinomially distributed random number vectors
  # from density function w1*multi1 + w2*multi2 + ... + wm*multim
  #
  # Args:
  #   n:    number of random vectors to draw.
  #   size: integer specifying the total number of objects that 
  #         are put into K boxes in the typical multinomial experiment. 
  #         where K is dim(p)[1]   
  #   prob:	a probability matrix of size K by m, specifying the 
  #         probability for the K classes for each of the m components 
  #         of the mixture
  #   w:    vector of weights for each of the multinomials used in the 
  #         mixture. If w == NULL, equal weights are used.
  # Returns:
  #   An integer K * n matrix where each column is a random vector
  #   generated according to the desired multionimial law and hence 
  #   summing to a size.
  
  p <- dim(prob)[2]
  wt <- if(is.null(w)) rep(1,p)/p else w
  c.out <- matrix(nrow=dim(prob)[1],ncol=n)
  for (i in 1:n) {
    v = matrix(data=0,nrow=dim(prob)[1],ncol=1)
    for (j in 1:size) {
      v = v + rmultinom(1,1,prob[,which(rmultinom(1,1,wt)==1)])
    }
    c.out[,i] <- v
  }
  return(c.out)
}


multinomResample=function(data,num=1){
  # Bootstrap resample of a data frame containing counts data
  # utilizing multinomial distributions.
  #
  # Args:
  #   data:   a data frame containing integer count data.
  #   num:    the number of bootstrap resamples to be generated.
  # Returns:
  #   A list of length equal to num containing data frames of
  #   bootstrap resamples of data.
  
  pr=as.matrix(data/sum(data))
  datlist=lapply(rep(1,num),RMultiMix,size=sum(data),prob=pr)
  return(datlist)
}


RAD=function(data,control,n=1000){
  # Performs RAD source selection for microbial source tracking for
  # rarefied abundance counts.
  #
  # Args:
  #   data:     array containing a single test/sink sample.
  #   control:  data frame containing all control/source samples.
  #   n:        the number of bootstrap resamples to generate for
  #             calculating boostrap estimates and test statistics.
  # Returns:
  #   A list of length two, the first element contains the RAD proportions
  #   for each control/source sample across all bootstrap permutations, the 
  #   second element contains a matrix of generated statistics and p-values
  #   Including the mean RAD score (RM) and its standard deviation (RSD), the
  #   mean RAD proportion (RP) and its standard deviation (RPSD), the lower
  #   bound of the 95% bootstrap interval (RLB) and the associated one-sided
  #   p-value (RPV), and an indicator for whether the RAD is significantly
  #   larger than zero (RI).
  
  require(robCompositions)
  dataset=control
  dataset$test=data
  dataset=ZeroFilter(dataset)
  nc=ncol(dataset)
  control=dataset[,-nc]
  #Generate New Samples
  listRS=as.list(rep(NA,n))
  for(i in 1:n){
    listRS[[i]]=control
  }
  print("bootstrap start")
  
  sampleRep=function(Data){
    # Internal function for apply integrated multinomResample.
    #
    # Args:
    #   Data:   a data frame containing integer count data.
    # Returns:
    #   A list of length one containing a data frame of
    #   one bootstrap resample of Data.
    
    remat=matrix(unlist(apply(Data,2,multinomResample)),nrow=nrow(Data))
    colnames(remat)=colnames(Data)
    return(remat)
  }
  
  listRS=lapply(listRS,sampleRep)
  print("bootstrap finished")
  
  
  #Create Data Sets
  combineDF2=function(data1,data2){
    # Internal function for apply integrated merging of test
    # and control samples.
    #
    # Args:
    #   data1:  data frame containing control samples
    #   data2:  array containing a single test sample
    # Returns:
    #   A data frame containing both control and test samples.
    
    comdata=data.frame(data1,"test"=data2)
    names(comdata)=c(colnames(data1),"test")
    return(comdata)
  }
  
  listData=lapply(listRS,combineDF2,data2=dataset$test)
  #Cluster and draw out cophenetic matrices
  listClust=listData
  
  for(i in 1:n){listClust[[i]]=listClust[[i]]+0.00001}
  aitchdist=function(Matrix){
    # Internal function for apply integrated aDist.
    #
    # Args:
    #   Matrix:   a matrix containing bootstrap resample dataset.
    # Returns:
    #   An Aitchison distance matrix for the given bootstrap resample.
    
    ADM=aDist(t(Matrix))
    return(ADM)
  }
  print("dist start")
  listClust=lapply(listClust,aitchdist)
  print("dist done")
  
  Outcome=as.list(rep(NA,length(methods)))
  
  intRAD=function(data){
    # Internal function for calculating RAD score using data
    # frame structure generated by combineDF2.
    #
    # Args:
    #   data:   A distance matrix for all control samples and
    #           a single test sample.
    # Returns:
    #   Array of RAD scores for all control samples.
    
    mmrd=function(array){
      # Internal function for calculating RAD score for a single
      # control sample.
      #
      # Args:
      #   array:  an array containing distances to all control
      #           samples and a single test sample.
      # Returns:
      #   RAD score for the given control sample.
      narray=array[which(array!=0)]
      nn=length(narray)
      return(mean(narray[-nn] - narray[nn])/mean(narray[-nn]))
    }
    
    data=as.matrix(data)
    n=nrow(data)
    smin=function(data){
      return(sort(data)[2])
    }
    smins=apply(data[-n,-n],2,smin)
    return(apply(data[,-n],2,mmrd))  
  }
  listRel=lapply(listClust,intRAD)
  
  #Calculate Estimates and Test Statistics from Bootstrap Resamples
  Rel=t(matrix(unlist(listRel),nrow=length(listRel[[1]])))
  RM=colMeans(Rel)
  RSD=apply(Rel,2,sd)
  RLB=RM-qnorm(1-0.05/10)*(RSD+0.001)/sqrt(n)
  z=RM/((RSD+0.001)/sqrt(n))
  RPV=pnorm(z)
  RI=(RLB>0)*1
  Draw=t(t(Rel)*RI)
  Drawsums=apply(Draw,1,function(x) sum(x))
  Draw=Draw/Drawsums
  RP=RM*(RI==1)
  RP=apply(Draw,2,mean)
  RPSD=apply(Draw,2,sd)
  RR=rbind(RM,RSD,RP,RPSD,RLB,RPV,RI)
  colnames(RR)=names(listRel[[1]])
  colnames(Draw)=names(listRel[[1]])
  
  #Combines data into a list
  Outcome=list(Draw,RR)
  return(Outcome)
}

