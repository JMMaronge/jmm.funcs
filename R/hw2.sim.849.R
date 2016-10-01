#' A function to do question 2 from homework 2 in STAT 849
#'
#' This function simulates data from  a distribution and performs a one sample T-test then outputs coverage probability and type 1 error rate 
#' @param h0.mean The mean under H0
#' @param num.samp The number of samples in each repetition     
#' @param num.reps The number of repetitions   
#' @param samp.dist The distribution to sample from, either "normal", "chisq", or "t"   
#' @param alpha Alpha level for test
#' @param times.scalar This is added for part g. If true, the function randomly samples a 5 or 10 for each repetition and multiples that by the first observation in that repetition's sample.
#' @param ... Additional options to specify for the sampling distribution
#' @keywords Classwork, STAT849, simulations
#' @export

hw2.sim.849<-function(h0.mean, num.samp,num.reps,samp.dist, times.scalar=FALSE, alpha=0.05,  ...){
  
  nsamp<-num.samp# number of samples in each rep
  nrep<-num.reps# number of repitions
  mu.cover<-vector(length = nrep)
  ci<-matrix(NA, nrow=nrep,ncol=2)
  samps<-matrix(NA,nrow=nrep,ncol=nsamp) #matrix of simulated data
  for(i in 1:nrep){
    if(samp.dist=="normal"){
      samps[i,]<- rnorm(nsamp, ...)} # generating from normal distribution
    else if(samp.dist=="t"){
      samps[i,]<-rt(nsamp, ...)} # generating from t distibution
    else if(samp.dist=="chisq"){ #generating from Chi-Squared
      samps[i,]<-rchisq(nsamp, ...)
    }
    else stop("This sampling distribution is not supported!") # check for samp.dist
    if(times.scalar==TRUE){ # for part g of the problem
      s<-c(5,10)
      pick<-sample(s,1)
      samps[i,1]<-samps[i,1]*pick}
    else{}
    test<-t.test(samps[i,], mu=h0.mean, conf.level = (1-alpha)) # performing t.test}
    mu.cover[i]<-(test$conf.int[1]<=h0.mean&h0.mean<=test$conf.int[2])
    ci[i,]<-test$conf.int
  }
  output<-vector(mode="list",length=3)
  names(output)<-c("cov.prob","type1.err.rate","CIs")
  output$cov.prob<- sum(mu.cover==TRUE)/nrep #calculating coverage probability
  output$type1.err.rate<- sum(mu.cover==FALSE)/nrep #calculates the type 1 error rate,
  #this works because if the 95% CI does not include 0, this implies p-val<.05
  output$CIs<-CIs
  return(output)
}