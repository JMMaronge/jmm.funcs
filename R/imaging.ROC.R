#' A function to create ROC plots from a list of prediction images for segmentation
#'
#' This function creates different types of ROC plots from prediction images for segmentation. It also calculates AUC and returns the data.frame used for the plot.
#' @param preds.list A list of prediction images as NIfTI objects for different subjects (can be of length 1)
#' @param mask.list A list of brain masks as NIfTI objects in the same order as preds.list
#' @param y.list A list of the gold standard segmentations as NIfTI objects in the same order as preds.list. 
#' @param vec.id A vector of ID's (as characters) corresponding to preds.list
#' @param type There are 3 options: "subject", "total", or "group". "subject" creates a seperate ROC for each subject, "total" creates one overall ROC, and by group makes a seperate ROC for each group as determined by the group.list argument.
#' @param group.list A list of sub-lists. This is only necessary if type="group". This should be a list of sublists where each sublist is the IDs of the subjects wanted in each group. For example, list(list("1","2"),list("3","4")) would create two ROCs, one with subjects "1" and "2" and the other with subjects "3" and "4".
#' @param max.fpr A scalar between (0,1). Use this to create partial ROCs instead of full ROCs. The ROC will be cut-off at the max false positve rate (max.fpr) specified.
#' @keywords ROC, AUC, NIfTI
#' @export
#' @examples
#' imaging.ROC()



imaging.ROC <-function(preds.list, mask.list, y.list,vec.id,type,group.list,max.fpr){
  require(ggplot2)
  require(data.table)
  if(type=="subject"){
    roc.subj.df<-vector(mode="list",length=length(preds.list))
    auc.subj.df<-vector(mode="list",length=length(preds.list))
    for(i in 1:length(preds.list)){
      print(paste0("Starting Subject ",i))
      preds<-preds.list[[i]]
      mask<-mask.list[[i]]
      y<-y.list[[i]]
      temp.dat<-data.frame(pred=preds[mask==1],les=y[mask==1])
      threshold <- seq(from = 1, to = 0, by= -1/1000)
      sens <-c()
      spec <-c()
      j<-1
      repeat{
        thresh <- threshold[j]
        tab <- matrix(c(sum(temp.dat$pred>=thresh & temp.dat$les==1),sum(temp.dat$pred<thresh & temp.dat$les==1),sum(temp.dat$pred>=thresh & temp.dat$les==0),sum(temp.dat$pred<thresh & temp.dat$les==0)),nrow=2,ncol=2)
        temp.sens <- (tab[1,1]/(tab[1,1]+tab[2,1]))
        temp.spec <- (tab[2,2]/(tab[2,2]+tab[1,2]))
        if(1-temp.spec>max.fpr){break}
        sens[j] <- temp.sens
        spec[j] <- temp.spec
        j<-j+1
      }
      # true.pos<- colSums(outer(temp.dat$pred,threshold,">=")&temp.dat$les==1)
      # false.neg<- colSums(outer(temp.dat$pred,threshold,"<")&temp.dat$les==1)
      # false.pos<- colSums(outer(temp.dat$pred,threshold,">=")&temp.dat$les==0)
      # true.neg<- colSums(outer(temp.dat$pred,threshold,"<")&temp.dat$les==0)
      # sens<-true.pos/(true.pos+false.neg)
      # spec<-true.neg/(false.pos+true.neg) 
      #Code works, but matrix too big for memory, alternative to for-loop above
      temp.ROC.data <- data.frame(ROCx = 1-spec, ROCy=sens,Threshold=threshold[1:(j-1)])  
      temp.ROC.data <- temp.ROC.data[order(temp.ROC.data$ROCx,temp.ROC.data$ROCy),]
      temp.ROC.data$Subject<-rep(vec.id[i],length(temp.ROC.data$ROCx))
      roc.subj.df[[i]]<-temp.ROC.data
      k <- 2:nrow(temp.ROC.data)
      temp.AUC <- data.frame(A=(temp.ROC.data$ROCx[k] - temp.ROC.data$ROCx[k - 1]) %*% (temp.ROC.data$ROCy[k] + temp.ROC.data$ROCy[k - 1])/2)
      temp.AUC$Subject<-vec.id[i]
      auc.subj.df[[i]]<-temp.AUC
    }
    ROC.data<-rbindlist(roc.subj.df)
    AUC.data<-rbindlist(auc.subj.df)
    ROC.plot <- ggplot(data=ROC.data,aes(x = ROCx, y=ROCy,colour=Subject)) +geom_point()+xlab("False Positive Rate (1-Specificity)")+ylab("True Positive Rate (Sensitivity)") + 
      ggtitle('ROC Curves \n for Each Subject') +
      theme(plot.title = element_text(lineheight=.8, face="bold"))+geom_abline (intercept = 0, slope = 1, colour="Red")
    roc.out<-vector(mode = "list",length = 3)
    roc.out$roc.data<-ROC.data
    roc.out$roc.plot<-ROC.plot
    roc.out$auc<-AUC.data
    return(roc.out) }  
  else if(type=="total"){
    subj.dat<-vector(mode="list",length = length(preds.list))
    for(i in 1:length(preds.list)){
      print(paste0("Starting Subject ",i))
      preds<-preds.list[[i]]
      mask<-mask.list[[i]]
      y<-y.list[[i]]
      subj.dat[[i]]<-data.frame(pred=preds[mask==1],les=y[mask==1])
    }
    total.dat<-rbindlist(subj.dat)
    print("Starting total calculations")
    threshold <- seq(from = 1, to = 0, by= -1/1000)
    sens <-c()
    spec <-c()
    j<-1
    repeat{
      thresh <- threshold[j]
      tab <- matrix(c(sum(total.dat$pred>=thresh & total.dat$les==1),sum(total.dat$pred<thresh & total.dat$les==1),sum(total.dat$pred>=thresh & total.dat$les==0),sum(total.dat$pred<thresh & total.dat$les==0)),nrow=2,ncol=2)
      temp.sens <- (tab[1,1]/(tab[1,1]+tab[2,1]))
      temp.spec <- (tab[2,2]/(tab[2,2]+tab[1,2]))
      if(1-temp.spec>max.fpr){break}
      sens[j] <- temp.sens
      spec[j] <- temp.spec
      j<-j+1
    }
    ROC.data <- data.frame(ROCx = 1-spec, ROCy=sens,Threshold=threshold[1:(j-1)])
    ROC.data <- ROC.data[order(ROC.data$ROCx,ROC.data$ROCy),]
    k <- 2:nrow(ROC.data)
    AUC.data <- data.frame(A=(ROC.data$ROCx[k] - ROC.data$ROCx[k - 1]) %*% (ROC.data$ROCy[k] + ROC.data$ROCy[k - 1])/2)
    ROC.plot <- ggplot(data=ROC.data,aes(x = ROCx, y=ROCy)) +geom_point()+xlab("False Positive Rate (1-Specificity)")+ylab("True Positive Rate (Sensitivity)") +
      ggtitle('ROC Curve') +
      theme(plot.title = element_text(lineheight=.8, face="bold"))+geom_abline (intercept = 0, slope = 1, colour="Red")
    roc.out<-vector(mode = "list",length = 3)
    roc.out$roc.data<-ROC.data
    roc.out$roc.plot<-ROC.plot
    roc.out$auc<-AUC.data
    return(roc.out)
  }
  else if(type=="group"){#### note list of subj id must be two sublists for by group
    subj.dat<-vector(mode="list",length = length(preds.list))
    for(i in 1:length(preds.list)){
      preds<-preds.list[[i]]
      mask<-mask.list[[i]]
      y<-y.list[[i]]
      temp.dat<-data.frame(pred=preds[mask==1],les=y[mask==1])
      temp.dat$Subject<-rep(vec.id[i],length(temp.dat$pred))
      subj.dat[[i]]<-temp.dat}
    total.dat<-rbindlist(subj.dat)
    
    
    group.ROC.list<-vector(mode = "list",length = length(group.list))
    group.AUC.list<-vector(mode = "list",length = length(group.list))
    for(i in 1:length(group.list)){
      print(paste0("Starting Group ",i))
      group<-group.list[[i]]
      group.dat<-total.dat[total.dat$Subject %in% group.list[[i]]]
      threshold <- seq(from = 1, to = 0, by= -1/1000)
      sens <-c()
      spec <-c()
      j<-1
      repeat{
        thresh <- threshold[j]
        tab <- matrix(c(sum(group.dat$pred>=thresh & group.dat$les==1),sum(group.dat$pred<thresh & group.dat$les==1),sum(group.dat$pred>=thresh & group.dat$les==0),sum(group.dat$pred<thresh & group.dat$les==0)),nrow=2,ncol=2)
        temp.sens <- (tab[1,1]/(tab[1,1]+tab[2,1]))
        temp.spec <- (tab[2,2]/(tab[2,2]+tab[1,2]))
        if(1-temp.spec>max.fpr){break}
        sens[j] <- temp.sens
        spec[j] <- temp.spec
        j<-j+1
      }
      temp.ROC.data <- data.frame(ROCx = 1-spec, ROCy=sens,Threshold=threshold[1:(j-1)])
      temp.ROC.data <- temp.ROC.data[order(temp.ROC.data$ROCx,temp.ROC.data$ROCy),]
      temp.ROC.data$Group<-rep(i,length(temp.ROC.data$ROCx))
      group.ROC.list[[i]]<-temp.ROC.data
      k <- 2:nrow(temp.ROC.data)
      temp.AUC <- data.frame(A=(temp.ROC.data$ROCx[k] - temp.ROC.data$ROCx[k - 1]) %*% (temp.ROC.data$ROCy[k] + temp.ROC.data$ROCy[k - 1])/2)
      temp.AUC$Group<-i
      group.AUC.list[[i]]<-temp.AUC
    }
    ROC.data<-rbindlist(group.ROC.list)
    AUC.data<-rbindlist(group.AUC.list)
    ROC.plot<-ggplot(data=ROC.data,aes(x = ROCx, y=ROCy,colour=factor(Group))) +geom_point()+xlab("False Positive Rate (1-Specificity)")+ylab("True Positive Rate (Sensitivity)") +
      ggtitle('ROC Curves \n for Each Group') +
      theme(plot.title = element_text(lineheight=.8, face="bold"))+geom_abline (intercept = 0, slope = 1, colour="Red")
    roc.out<-vector(mode = "list",length = 3)
    roc.out$roc.data<-ROC.data
    roc.out$roc.plot<-ROC.plot
    roc.out$auc<-AUC.data
    roc.out$groups<- group.list
    return(roc.out)
  }
  else{print("ERROR: argument for 'type' not recognized!")
  }
}