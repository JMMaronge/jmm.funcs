#' A function to do question 2 from homework 2 in STAT 849
#'
#' This function simulates data from  a distribution and performs a one sample T-test then outputs coverage probability and type 1 error rate 
#' @param h0.mean The mean under H0
#' @param num.samp The number of samples in each repetition     
#' @param num.reps The number of repetitions   
#' @param samp.dist The distribution to sample from, either "normal", "chisq", or "t"   
#' @param times.scalar This is added for part g. If true, the function randomly samples a 5 or 10 for each repetition and multiples that by the first observation in that repetition's sample.
#' @param ... Additional options to specify for the sampling distribution
#' @keywords Classwork, STAT849, simulations
#' @export

hw2.sim.849<-function(h0.mean, num.samp,num.reps,samp.dist, times.scalar=FALSE,  ...){
  
  nsamp<-num.samp# number of samples in each rep
  nrep<-num.reps# number of repitions
  mu.cover<-vector(length = nrep)
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
    test<-t.test(samps[i,], mu=h0.mean) # performing t.test}
    mu.cover[i]<-prod((test$conf.int - h0.mean))<0 # if this is true, 0 is within the ci, o.w. false
  }
  output<-vector(mode="list",length=2)
  names(output)<-c("cov.prob","type1.err.rate")
  output$cov.prob<- sum(mu.cover==TRUE)/nrep #calculating coverage probability
  output$type1.err.rate<- sum(mu.cover==FALSE)/nrep #calculates the type 1 error rate,
  #this works because if the 95% CI does not include 0, this implies p-val<.05
  return(output)
}