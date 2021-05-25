# rm(list=ls())
Cross.Validation<-function(Spec     			  ## Spektren
                           ,Labels				  ## Labels-Matrix bzw. Regressantenvektor
                           ,Batch=NULL			## Batchvektor fuer Batchwise-Validation
                           ,Batch.N=8				## Number for Batch-out-validation (ignored if Batch=NULL)
                           ,N=8					  ## Number for Crossvalidation (ignored if Batch supplied)
                           ,Equalize.Batch=FALSE		## Equalize Batch nur für Batchvalidation & Classification (Currently Ignored)
                           ,Print.Output=TRUE
                           ,...)
{
  if (is.null(Batch))      # there's no batch vector based on which for CV
  {
    if (Print.Output)
    {
      print(paste("CV",N))
    }
    FOLD<-vector("list",N)
    SAMPLE<-sample(1:NROW(Spec),NROW(Spec))
    REST<-NROW(Spec)%%N       # rest after partition
    PART<-NROW(Spec)%/%N      # partition
    for (i in 1:N)            # distribute the rest spectra
    {
      FOLD[[i]]<-SAMPLE[1:PART+(i-1)*PART]
      if (REST>0)
      {
        FOLD[[i]]<-c(FOLD[[i]],SAMPLE[N*PART+REST])
        REST<-REST-1
      }
    }
  }
  else	
  {
    LEV<-unique(Batch)
    if (Batch.N>1)
    {
      if (Print.Output)
      {
        print(paste("BV",N,"by",Batch.N))
      }
      LEVELS<-vector("list",Batch.N)
      Index<-1:length(LEV)
      for (i in 1:Batch.N)
      {
        LEVELS[[i]]<-LEV[which(Index%%Batch.N==(i-1))]  # to partion "batch" in Batch.N parts
      }
      FOLD<-vector("list",Batch.N)
      for (i in 1:Batch.N)
      {
        FOLD[[i]]<-which(is.element(el=Batch,set=LEVELS[[i]]))
      }
    }
    else
    {
      N<-length(LEV)
      if (Print.Output)
      {
        print(paste("BV",N))
      }
      FOLD<-vector("list",N)
      for (i in 1:N)
      {
        FOLD[[i]]<-which(Batch==LEV[i])
      }
    }
  }
  
  return(FOLD)
}

Make.classify.pca <- function(Spec, Labels, Batch=NULL, Test, testLabels, Class.fn=svm, Pred.fn=predict, ncomp=3:50, ...)
{
  FOLD <- Cross.Validation(Spec=Spec, Labels=Labels, Batch=Batch, Print.Output=FALSE)
  
  Labels <- factor(Labels, levels=unique(Labels))
  
  N <- length(FOLD)
  SENS <- c()
  for (i in 1:N)
  {
    mSens <- c()
    INDEX<-1:NROW(Spec)
    INDEX<-INDEX[-FOLD[[i]]]
    PCA.Model <- prcomp(Spec[INDEX,], retx = TRUE, center = FALSE, scale. = FALSE)
    pscores <- predict(PCA.Model, Spec[FOLD[[i]],,drop=F])
    for(nPC in ncomp)
    {
      tscores <- PCA.Model$x[,1:nPC]
      if(!is.matrix(tscores)) tscores <- as.matrix(tscores)
      if(!is.matrix(pscores)) pscores <- as.matrix(pscores)
      
      CLASS<-Class.fn(tscores, Labels[INDEX], ...)
      TMP.p<-Pred.fn(CLASS, pscores[,1:nPC])
      
      True <- as.character(Labels)[FOLD[[i]]]    # true values
      Pred <- as.character(TMP.p)
      
      True <- factor(True, levels=levels(Labels))
      Pred <- factor(Pred, levels=levels(Labels))
      
      SUMM <- confusionMatrix(Pred, True)
      if(length(unique(Labels))>2)
        mSens <- c(mSens, mean(SUMM$byClass[, 'Sensitivity'], na.rm=TRUE))
      else
        mSens <- c(mSens, mean(SUMM$byClass[c('Sensitivity', 'Specificity')]))
    }
    colnames(SENS) <- NULL
    SENS <- cbind(SENS, mSens)
  }
  
  optPC <- ncomp[which.max(rowMeans(SENS, na.rm=TRUE))]
  
  PCA.Model <- prcomp(Spec, retx = TRUE, center = FALSE, scale. = FALSE)
  
  tscores <- PCA.Model$x[,1:optPC]
  if(!is.matrix(tscores)) tscores <- as.matrix(tscores)
  
  pscores <- predict(PCA.Model, Test)[,1:optPC]
  if(!is.matrix(pscores)) pscores <- as.matrix(pscores)
  
  CLASS <- Class.fn(tscores, Labels,...)
  TMP.p <- Pred.fn(CLASS, pscores)
  
  True <- as.character(testLabels)    # true values
  Pred <- as.character(TMP.p)
  
  True <- factor(True, levels=levels(Labels))
  Pred <- factor(Pred, levels=levels(Labels))
  
  SUMM <- confusionMatrix(Pred, True)
  if(length(unique(Labels))>2)
    mSens <- mean(SUMM$byClass[, 'Sensitivity'], na.rm=TRUE)
  else
    mSens <- mean(SUMM$byClass[c('Sensitivity', 'Specificity')])
  
  return(list(SENS=SENS, testSens=mSens, optPC=optPC))
}

Make.classify.pca1 <- function(Spec, Labels, Batch=NULL, Test, testLabels, Class.fn=svm, Pred.fn=predict, ncomp=3:50, ...)
{
  FOLD <- Cross.Validation(Spec=Spec, Labels=Labels, Batch=Batch, Print.Output=FALSE)
  
  Labels <- factor(Labels, levels=unique(Labels))
  
  N<-length(FOLD)
  
  PCA.Model <- prcomp(Spec, retx = TRUE, center = FALSE, scale. = FALSE)
  tscores <- PCA.Model$x
  if(!is.matrix(tscores)) tscores <- as.matrix(tscores)
  
  SENS <- c()
  for (i in 1:N)
  {
    mSens <- c()
    INDEX<-1:NROW(Spec)
    INDEX<-INDEX[-FOLD[[i]]]
    for(nPC in ncomp)
    {
      CLASS<-Class.fn(as.matrix(tscores[INDEX,1:nPC]), Labels[INDEX],...)
      TMP.p <- Pred.fn(CLASS, as.matrix(tscores[FOLD[[i]],1:nPC]))
      True<-as.character(Labels)[FOLD[[i]]]    # true values
      Pred<-as.character(TMP.p)
      
      True <- factor(True, levels=levels(Labels))
      Pred <- factor(Pred, levels=levels(Labels))
      
      SUMM <- confusionMatrix(Pred, True)
      if(length(unique(Labels))>2)
        mSens <- c(mSens, mean(SUMM$byClass[, 'Sensitivity'], na.rm=TRUE))
      else
        mSens <- c(mSens, mean(SUMM$byClass[c('Sensitivity', 'Specificity')]))
    }
    colnames(SENS) <- NULL
    SENS <- cbind(SENS, mSens)
  }
  optPC <- ncomp[which.max(rowMeans(SENS, na.rm=TRUE))]
  
  tscores <- PCA.Model$x[,1:optPC]
  if(!is.matrix(tscores)) tscores <- as.matrix(tscores)
  
  pscores <- predict(PCA.Model, Test)[,1:optPC]
  if(!is.matrix(pscores)) pscores <- as.matrix(pscores)
  
  CLASS<-Class.fn(tscores, Labels,...)
  TMP.p <-Pred.fn(CLASS, pscores)
  
  True<-as.character(testLabels)    # true values
  Pred<-as.character(TMP.p)
  
  
  True <- factor(True, levels=levels(Labels))
  Pred <- factor(Pred, levels=levels(Labels))
  
  SUMM <- confusionMatrix(Pred, True)
  if(length(unique(Labels))>2)
    mSens <- mean(SUMM$byClass[, 'Sensitivity'], na.rm=TRUE)
  else
    mSens <- mean(SUMM$byClass[c('Sensitivity', 'Specificity')])
  
  
  return(list(SENS=SENS, testSens=mSens, optPC=optPC))
}

Make.classify.pls <- function(Spec, Labels, Batch=NULL, Test, testLabels, Class.fn=svm, Pred.fn=predict, ncomp=3:50, ...)
{
  require(pls)
  FOLD <- Cross.Validation(Spec=Spec, Labels=Labels, Batch=Batch, Print.Output=FALSE)
  
  N<-length(FOLD)
  Y.matrix <- matrix(0, nrow=nrow(Spec), ncol=length(unique(Labels)))
  for(j in 1:ncol(Y.matrix))
  {
    ix <- which(Labels==unique(Labels)[j])
    Y.matrix[ix,j] <- 1
  }
  
  Labels <- factor(Labels, levels=unique(Labels))
  
  SENS <- c()
  for (i in 1:N)
  {
    mSens <- c()
    INDEX<-1:NROW(Spec)
    INDEX<-INDEX[-FOLD[[i]]]
    
    data <- list(Response=Y.matrix[INDEX], Variable=Spec[INDEX,])
    
    for(nCp in ncomp)
    {
      Model <- plsr(Response~Variable, ncomp=nCp, data=data, scale=TRUE, validation="none", y=T)
      tscores <- scores(Model)
      if(!is.matrix(tscores)) tscores <- as.matrix(tscores)
      scores.m <- c()
      for(j in 1:nrow(tscores))
      {
        scores.m <- rbind(scores.m, tscores[j,])
      }
      tscores <- as.matrix(scores.m)
      
      CLASS<-Class.fn(x=tscores, Labels[INDEX],...)
      
      pscores <- predict(Model, Spec[FOLD[[i]],,drop=F], type="scores")
      if(!is.matrix(pscores)) pscores <- as.matrix(pscores)
      scores.m <- c()
      for(j in 1:nrow(pscores))
      {
        scores.m <- rbind(scores.m, pscores[j,])
      }
      pscores <- as.matrix(scores.m)
      TMP.p <- Pred.fn(CLASS, pscores)
      
      True<-as.character(Labels)[FOLD[[i]]]    # true values
      Pred<-as.character(TMP.p)
      
      True <- factor(True, levels=levels(Labels))
      Pred <- factor(Pred, levels=levels(Labels))
      
      SUMM <- confusionMatrix(Pred, True)
      if(length(unique(Labels))>2)
        mSens <- c(mSens, mean(SUMM$byClass[, 'Sensitivity'], na.rm=TRUE))
      else
        mSens <- c(mSens, mean(SUMM$byClass[c('Sensitivity', 'Specificity')]))
    }
    colnames(SENS) <- NULL
    SENS <- cbind(SENS, mSens)
  }
  optCp <- ncomp[which.max(rowMeans(SENS, na.rm=TRUE))]
  
  data <- list(Response=Y.matrix, Variable=Spec)
  Model <- plsr(Response~Variable, ncomp=optCp, data=data, scale=TRUE, validation="none", y=T)
  tscores <- scores(Model)
  if(!is.matrix(tscores)) tscores <- as.matrix(tscores)
  scores.m <- c()
  for(j in 1:nrow(tscores))
  {
    scores.m <- rbind(scores.m, tscores[j,])
  }
  tscores <- as.matrix(scores.m)
  
  CLASS <- Class.fn(x=tscores, Labels,...)
  
  pscores <- predict(Model, Test, type="scores")
  if(!is.matrix(pscores)) pscores <- as.matrix(pscores)
  scores.m <- c()
  for(j in 1:nrow(pscores))
  {
    scores.m <- rbind(scores.m, pscores[j,])
  }
  pscores <- as.matrix(scores.m)
  TMP.p <- Pred.fn(CLASS, pscores)
  
  True<-as.character(testLabels)    # true values
  Pred<-as.character(TMP.p)
  
  True <- factor(True, levels=levels(Labels))
  Pred <- factor(Pred, levels=levels(Labels))
  
  SUMM <- confusionMatrix(Pred, True)
  if(length(unique(Labels))>2)
    mSens <- mean(SUMM$byClass[, 'Sensitivity'], na.rm=TRUE)
  else
    mSens <- mean(SUMM$byClass[c('Sensitivity', 'Specificity')])
  
  
  return(list(SENS=SENS, testSens=mSens, optCp=optCp))
}

Make.classify.pls1 <- function(Spec, Labels, Batch=NULL, Test, testLabels, Class.fn=svm, Pred.fn=predict, ncomp=3:50, ...)
{
  require(pls)
  FOLD <- Cross.Validation(Spec=Spec, Labels=Labels, Batch=Batch, Print.Output=FALSE)
  
  N<-length(FOLD)
  
  Y.matrix <- matrix(0, nrow=nrow(Spec), ncol=length(unique(Labels)))
  for(j in 1:ncol(Y.matrix))
  {
    ix <- which(Labels==unique(Labels)[j])
    Y.matrix[ix,j] <- 1
  }
  data <- list(Response=Y.matrix, Variable=Spec)
  
  Labels <- factor(Labels, levels=unique(Labels))
  
  SENS <- c()
  for(nCp in ncomp)
  {
    Model <- plsr(Response~Variable, ncomp=nCp, data=data, scale=TRUE, validation="none", y=T)
    tscores <- scores(Model)
    if(!is.matrix(tscores)) tscores <- as.matrix(tscores)
    scores.m <- c()
    for(j in 1:nrow(tscores))
    {
      scores.m <- rbind(scores.m, tscores[j,])
    }
    tscores <- as.matrix(scores.m)
    mSens <- c()
    for (i in 1:N)
    {
      INDEX<-1:NROW(Spec)
      INDEX<-INDEX[-FOLD[[i]]]   
      
      CLASS <- Class.fn(as.matrix(tscores[INDEX,]), Labels[INDEX],...)
      
      TMP.p <- Pred.fn(CLASS, as.matrix(tscores[FOLD[[i]],]))
      
      True<-as.character(Labels)[FOLD[[i]]]    # true values
      Pred<-as.character(TMP.p)
      
      True <- factor(True, levels=levels(Labels))
      Pred <- factor(Pred, levels=levels(Labels))
      
      SUMM <- confusionMatrix(Pred, True)
      if(length(unique(Labels))>2)
        mSens <- c(mSens, mean(SUMM$byClass[, 'Sensitivity'], na.rm=TRUE))
      else
        mSens <- c(mSens, mean(SUMM$byClass[c('Sensitivity', 'Specificity')]))
    }
    SENS <- rbind(SENS, mSens)
  }
  
  optCp <- ncomp[which.max(rowMeans(SENS, na.rm=TRUE))]
  
  Model <- plsr(Response~Variable, ncomp=optCp, data=data, scale=TRUE, validation="none", y=T)
  tscores <- scores(Model)
  if(!is.matrix(tscores)) tscores <- as.matrix(tscores)
  scores.m <- c()
  for(j in 1:nrow(tscores))
  {
    scores.m <- rbind(scores.m, tscores[j,])
  }
  tscores <- as.matrix(scores.m)
  
  CLASS <- Class.fn(x=tscores, Labels,...)
  
  pscores <- predict(Model, Test, type="scores")
  if(!is.matrix(pscores)) pscores <- as.matrix(pscores)
  scores.m <- c()
  for(j in 1:nrow(pscores))
  {
    scores.m <- rbind(scores.m, pscores[j,])
  }
  pscores <- as.matrix(scores.m)
  TMP.p <- Pred.fn(CLASS, pscores)
  
  True<-as.character(testLabels)    # true values
  Pred<-as.character(TMP.p)
  
  True <- factor(True, levels=levels(Labels))
  Pred <- factor(Pred, levels=levels(Labels))
  
  SUMM <- confusionMatrix(Pred, True)
  if(length(unique(Labels))>2)
    mSens <- mean(SUMM$byClass[, 'Sensitivity'], na.rm=TRUE)
  else
    mSens <- mean(SUMM$byClass[c('Sensitivity', 'Specificity')])
  
  
  return(list(SENS=SENS, testSens=mSens, optCp=optCp))
}

predLDA <- function(model, data)
{
  return(predict(model, data)$class)
}

predSVM <- function(model, data)
{
  return(predict(model, data))
}

predRF <- function(model, data)
{
  return(predict(model, data))
}

Plot.Accuracy <- function(Index, DATA, Names,YLIM=NULL,TYPE="b",MAIN="",...)
{
  
  if (is.null(YLIM))
  {
    YLIM<-range(c(DATA[[1]], DATA[[2]], DATA[[3]], DATA[[4]],
                  DATA[[5]], DATA[[6]], DATA[[7]], DATA[[8]]))
  }
  
  split.screen(figs=rbind(c(0,0.43,0,1),c(0.43,0.8,0,1),c(0.8,1,0,1)))
  screen(1)
  par(cex=1, mar=c(4.1, 4.1, 2.61, 0.51), cex.main=0.9, cex.axis=0.8, font.main=1)
  plot(Index,DATA[[1]],ylim=YLIM,type=TYPE,ylab="Mean Sensitivity"
       ,main="Inside Loop",xlab="Dimension Index", cex.lab=0.8, ...)
  points(Index,DATA[[2]],pch=4,type=TYPE)
  points(Index,DATA[[3]],col=2,type=TYPE)
  points(Index,DATA[[4]],pch=4,col=2,type=TYPE)
  screen(2)
  par(mar=c(4.1, 1.51, 2.61, 1.1),yaxt="n", cex.lab=1)
  plot(Index,DATA[[5]],ylim=YLIM,col=3,type=TYPE,ylab=""
       ,main="Outside Loop",xlab="Dimension Index",cex.lab=0.8)
  points(Index,DATA[[6]],pch=4,col=3,type=TYPE)
  points(Index,DATA[[7]],col=4,type=TYPE)
  points(Index,DATA[[8]],pch=4,col=4,type=TYPE)
  screen(3)
  par(xaxt="n",yaxt="n",bty="n",mar=rep(0,4))
  plot(0:10/10,0:10/10,type="n",xlab="",ylab="")
  mtext(side=3,line=-2,MAIN,font=2, cex=0.8)
  legend(x=0,y=0.8,legend=paste(c("PCA k-Replicate","PCA k-Fold","PLS k-Replicate","PLS k-Fold","PCA k-Replicate","PCA k-Fold","PLS k-Replicate","PLS k-Fold"),
                                sep=" ")
         ,pch=c(1,4,1,4,1,4,1,4),col=c(1,1,2,2,3,3,4,4), cex=0.7)
  close.screen(all.sc=T)
}

plot.accSd <- function(Index, DATA, SD, Names,YLIM=NULL,TYPE="l",MAIN="",...)
{
  
  if (is.null(YLIM))
  {
    YLIM<-range(c(DATA[[1]], DATA[[2]], DATA[[3]], DATA[[4]],
                  DATA[[5]], DATA[[6]], DATA[[7]], DATA[[8]]))
  }
  
  split.screen(figs=rbind(c(0,0.42,0,1),c(0.42,0.78,0,1),c(0.78,1,0,1)))
  screen(1)
  par(cex=1, mar=c(4.1, 4.1, 2.61, 0.51), cex.main=0.9, cex.axis=0.8, font.main=1)
  plot(Index,DATA[[1]],ylim=YLIM,type=TYPE,ylab="Mean Sensitivity"
       ,main="Inside-CV",xlab=expression(paste(italic("nPC"), ' (', italic("nLV"), ") ", sep='')), cex.lab=0.8, ...)
  polygon(c(Index, rev(Index)), c(DATA[[1]]+SD[[1]]/2, rev(DATA[[1]]-SD[[1]]/2)), col=rgb(0.7,0.7,0.7,0.5), border=FALSE)
  lines(Index, DATA[[1]])
  polygon(c(Index, rev(Index)), c(DATA[[2]]+SD[[2]]/2, rev(DATA[[2]]-SD[[2]]/2)), col=rgb(0.7,0.7,0.7,0.5), border=FALSE)
  lines(Index,DATA[[2]],lty=2)
  polygon(c(Index, rev(Index)), c(DATA[[3]]+SD[[3]]/2, rev(DATA[[3]]-SD[[3]]/2)), col=rgb(0.9,0.5,0.5,0.5), border=FALSE)
  lines(Index,DATA[[3]],lty=1,col=rgb(1,0,0))
  polygon(c(Index, rev(Index)), c(DATA[[4]]+SD[[4]]/2, rev(DATA[[4]]-SD[[4]]/2)), col=rgb(0.9,0.5,0.5,0.5), border=FALSE)
  lines(Index,DATA[[4]],lty=2,col=rgb(1,0,0))
  
  screen(2)
  par(mar=c(4.1, 1.51, 2.61, 1.1),yaxt="n", cex.lab=1)
  plot(Index,DATA[[5]],ylim=YLIM,col=rgb(0,1,0),type=TYPE,ylab=""
       ,main="Outside-CV",xlab=expression(paste(italic("nPC"), ' (', italic("nLV"), ") ", sep='')),cex.lab=0.8)
  polygon(c(Index, rev(Index)), c(DATA[[5]]+SD[[5]]/2, rev(DATA[[5]]-SD[[5]]/2)), col=rgb(0.5,0.9,0.5,0.5), border=FALSE)
  lines(Index,DATA[[5]],lty=1,col=rgb(0,1,0))
  polygon(c(Index, rev(Index)), c(DATA[[6]]+SD[[6]]/2, rev(DATA[[6]]-SD[[6]]/2)), col=rgb(0.5,0.9,0.5,0.5), border=FALSE)
  lines(Index,DATA[[6]],lty=2,col=rgb(0,1,0))
  polygon(c(Index, rev(Index)), c(DATA[[7]]+SD[[7]]/2, rev(DATA[[7]]-SD[[7]]/2)), col=rgb(0.5,0.5,0.9,0.5), border=FALSE)
  lines(Index,DATA[[7]],lty=1,col=rgb(0,0,1))
  polygon(c(Index, rev(Index)), c(DATA[[8]]+SD[[8]]/2, rev(DATA[[8]]-SD[[8]]/2)), col=rgb(0.5,0.5,0.9,0.5), border=FALSE)
  lines(Index,DATA[[8]],lty=2,col=rgb(0,0,1))
  
  screen(3)
  par(xaxt="n",yaxt="n",bty="n",mar=rep(0,4))
  plot(0:10/10,0:10/10,type="n",xlab="",ylab="")
  mtext(side=3,line=-2,MAIN,font=2, cex=0.8)
  legend(x=0,y=0.8,legend=paste(c("PCA k-Replicate","PCA k-Fold","PLS k-Replicate","PLS k-Fold","PCA k-Replicate","PCA k-Fold","PLS k-Replicate","PLS k-Fold"),
                                sep=" ")
         ,lty=c(1,2,1,2,1,2,1,2),col=c(1,1,2,2,3,3,4,4), cex=0.7)
  close.screen(all.sc=T)
}

plot.Box <- function(DATA, MAIN="",xtext=TRUE, method, mar=c(0,0,0,0), ...)
{
  par(mar=mar)
  clsLevel <- c("pca.bv.I", "pca.fv.I", "pca.bv.O", "pca.fv.O",
                "pls.bv.I", "pls.fv.I", "pls.bv.O", "pls.fv.O")
  data <- c()
  group <- c()
  for(i in 1:length(DATA))
  {
    data <- c(data, 100*DATA[[i]])
    group <- c(group, rep(clsLevel[i], length(DATA[i][[1]])))
  }
  data.plot <- cbind(factor(group, levels=clsLevel), data)
  colnames(data.plot) <- c("labels", "X")
  boxplot(formula=X~labels, data=data.plot, xlab='', ylim=c(0,110),
          xlab='', ylab="", names=NULL, xaxt='n',boxwex=0.2, main=MAIN, ...)
  if(length(MAIN)>0) mtext('Mean Sensitivity / %', side=2, line=2.5)
  for(i in 1:length(DATA))
  {
    lines(c(i,i),par('usr')[3:4], type="l", lty=2, col=rgb(0.8,0.8,0.8))
  }
  text(length(clsLevel)+0.3,par('usr')[4]-5, method)
  clsLevel1 <- c("pca.R.I", "pca.F.I", "pca.R.O", "pca.F.O",
                 "pls.R.I", "pls.F.I", "pls.R.O", "pls.F.O")
  if(xtext)
    text(1:length(clsLevel), y=par('usr')[3]-5, labels=clsLevel1, srt=-30, adj=0, xpd=NA, cex=1)
  box()
}

statTest <- function(DATA=list(PCA.LDA.BV.mb.I[, 97],
                               PCA.LDA.CV.mb.I[,97],
                               PCA.LDA.BV.mb.O[,97],
                               PCA.LDA.CV.mb.O[,97],
                               PLS.LDA.BV.mb.I[,97],
                               PLS.LDA.CV.mb.I[,97],
                               PLS.LDA.BV.mb.O[,97],
                               PLS.LDA.CV.mb.O[,97]),
                     item1=c(1,3,5,7), item2=c(2,4,6,8),
                     fHeight=3, colLine=rgb(1,0.5,0.5), colStar=rgb(1,0,0))
{
  for(i in 1:length(item1))
  {
    sigResult <- wilcox.test(x=100*DATA[item1[i]][[1]], y=100*DATA[item2[i]][[1]], alternative="less", conf.level = 0.95, exact=TRUE)
    pTest <- sigResult$p.value
    lines(c(item1[i], item2[i]), c(par('usr')[3]+fHeight+fHeight*(item1[i]%%2), par('usr')[3]+fHeight+fHeight*(item1[i]%%2)), col=colLine, lty=2)
    lines(c(item1[i],item1[i]), c(par('usr')[3], par('usr')[3]+fHeight+fHeight*(item1[i]%%2)), col=colLine, lty=2)
    lines(c(item2[i],item2[i]), c(par('usr')[3], par('usr')[3]+fHeight+fHeight*(item1[i]%%2)), col=colLine, lty=2)
    text(c(item1[i]+item2[i])*0.5, par('usr')[3]+fHeight+fHeight*(item1[i]%%2), labels=round(pTest,2), pos=3, offset=0.1,col=colStar)
  }
}