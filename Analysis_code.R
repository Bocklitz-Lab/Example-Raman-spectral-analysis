rm(list=ls())
while(dev.cur()!=1)
{
  dev.off()
}

forPath <- function(x){ x<-1 }
setwd(paste(dirname(dirname(getSrcDirectory(forPath))), '/', sep=''))

require(e1071)
require(MASS)
require(randomForest)
require(caret)

source(paste(getwd(), '/Scripts/functions_1.R', sep=''))
load(paste(getwd(), "/Data/DATA_dp_wc_bc.RData", sep=''))

# v_csv <- cbind(DATA.c$NAME.m, DATA.c$DATA.m, DATA.c$Bewertung.m[,1], DATA.c$Tissue.m)
# colnames(v_csv) <- c('Name', colnames(DATA.c$DATA.m), 'Annotation', 'Tissue')
# write.csv(v_csv, file=paste(getwd(), "/Data/DATA_dp_wc_bc.csv", sep=''))

ix <- which(DATA.c$Bewertung.m[, 1]=='Hp')
spec <- DATA.c$DATA.m[-ix, ]
labels <- DATA.c$Bewertung.m[-ix, 1]
ids <- DATA.c$NAME.m[-ix]

batches <- c()
for(i in 1:length(ids))
{
  batches <- c(batches, strsplit(ids[i], "_")[[1]][1])
}
tissues <- DATA.c$Tissue.m[-ix]
rm(DATA.c)

labels[which(labels%in%c("Normal"))] <- 'Normal'
labels[which(labels%in%c("Adenom","Karzinom"))] <- 'Abnormal'

vLevels <- c("Normal", "Abnormal")
### get wn region of interest
wn <- as.numeric(colnames(spec))
ixSel <- which(wn>600 & wn<1800)
wn <- wn[ixSel]
spec <- spec[, ixSel]/sqrt(rowSums(spec[, ixSel]^2))

### prepare data for model training
ixPrep <- which(tissues=="preparation")
spec <- spec[ixPrep,]
labels <- labels[ixPrep]
batches <- batches[ixPrep]

fCV <- Cross.Validation(Spec=specTrain, Labels=labels,
                        Batch=batches,
                        Print.Output=FALSE,
                        Batch.N=9)


# PCA+BV/CV Inside
PCA.LDA.BV.I <- c()
PCA.LDA.CV.I <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pca(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=lda, Pred.fn=predLDA)
  PCA.LDA.BV.I <- rbind(PCA.LDA.BV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))

  TMP <- Make.classify.pca(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=lda, Pred.fn=predLDA)
  PCA.LDA.CV.I <- rbind(PCA.LDA.CV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
}

# PLS+BV/CV Inside
PLS.LDA.BV.I <- c()
PLS.LDA.CV.I <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pls(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=lda, Pred.fn=predLDA)
  PLS.LDA.BV.I <- rbind(PLS.LDA.BV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))

  TMP <- Make.classify.pls(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=lda, Pred.fn=predLDA)
  PLS.LDA.CV.I <- rbind(PLS.LDA.CV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
}

# PCA+BV/CV Outside
PCA.LDA.BV.O <- c()
PCA.LDA.CV.O <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pca1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=lda,  Pred.fn=predLDA)
  PCA.LDA.BV.O <- rbind(PCA.LDA.BV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))

  TMP <- Make.classify.pca1(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=lda, Pred.fn=predLDA)
  PCA.LDA.CV.O <- rbind(PCA.LDA.CV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
}

# PLS+BV/CV Outside
PLS.LDA.BV.O <- c()
PLS.LDA.CV.O <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pls1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=lda, Pred.fn=predLDA)
  PLS.LDA.BV.O <- rbind(PLS.LDA.BV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))

  TMP <- Make.classify.pls1(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=lda, Pred.fn=predLDA)
  PLS.LDA.CV.O <- rbind(PLS.LDA.CV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
}

save(file=paste(getwd(), "/Results/LDA_Mice.RData", sep=''), list=c("PCA.LDA.BV.I",
                                                                         "PCA.LDA.CV.I",
                                                                         "PLS.LDA.BV.I",
                                                                         "PLS.LDA.CV.I",
                                                                         "PCA.LDA.BV.O",
                                                                         "PCA.LDA.CV.O",
                                                                         "PLS.LDA.BV.O",
                                                                         "PLS.LDA.CV.O"))

load(file=paste(getwd(), "/Results/LDA_Mice.RData", sep=''))
k=1
pdf('Mice_LDA.pdf', width=7, height=3)
plot.accSd(Index=3:50, DATA=list(PCA.LDA.BV.I[k, 1:48],
                                 PCA.LDA.CV.I[k, 1:48],
                                 PLS.LDA.BV.I[k, 1:48],
                                 PLS.LDA.CV.I[k, 1:48],
                                 PCA.LDA.BV.O[k, 1:48],
                                 PCA.LDA.CV.O[k, 1:48],
                                 PLS.LDA.BV.O[k, 1:48],
                                 PLS.LDA.CV.O[k, 1:48]), 
           SD=list(PCA.LDA.BV.I[k, 49:96],
                   PCA.LDA.CV.I[k, 49:96],
                   PLS.LDA.BV.I[k, 49:96],
                   PLS.LDA.CV.I[k, 49:96],
                   PCA.LDA.BV.O[k, 49:96],
                   PCA.LDA.CV.O[k, 49:96],
                   PLS.LDA.BV.O[k, 49:96],
                   PLS.LDA.CV.O[k, 49:96]), 
           Names,YLIM=NULL,TYPE="l",MAIN="Mice Data (LDA)")
dev.off()

pdf('box_Mice_LDA.pdf', width=5, height=3)
# nf <- layout(matrix(c(1,2),nrow=1,byrow=FALSE), widths=c(4,1))
plot.Box(DATA=list(PCA.LDA.BV.I[, 97],
                   PCA.LDA.CV.I[,97],
                   PCA.LDA.BV.O[,97],
                   PCA.LDA.CV.O[,97],
                   PLS.LDA.BV.I[,97],
                   PLS.LDA.CV.I[,97],
                   PLS.LDA.BV.O[,97],
                   PLS.LDA.CV.O[,97]), 
         MAIN="Mice Data (LDA)", xtext=TRUE, method="", mar=c(3,4,2,2), cex.main=0.8)

plot.Box(DATA=list(apply(PCA.LDA.BV.I[,1:48],1,max),
                   apply(PCA.LDA.CV.I[,1:48],1,max),
                   apply(PCA.LDA.BV.O[,1:48],1,max),
                   apply(PCA.LDA.CV.O[,1:48],1,max),
                   apply(PLS.LDA.BV.I[,1:48],1,max),
                   apply(PLS.LDA.CV.I[,1:48],1,max),
                   apply(PLS.LDA.BV.O[,1:48],1,max),
                   apply(PLS.LDA.CV.O[,1:48],1,max)),
         MAIN="", xtext=FALSE, method="", mar=c(3,4,2,2), at=(1:8)+0.3, add=TRUE,border=rgb(0.4,0.4,0.4))
# statTest(list(PCA.LDA.BV.I[, 97],
#               PCA.LDA.CV.I[,97],
#               PCA.LDA.BV.O[,97],
#               PCA.LDA.CV.O[,97],
#               PLS.LDA.BV.I[,97],
#               PLS.LDA.CV.I[,97],
#               PLS.LDA.BV.O[,97],
#               PLS.LDA.CV.O[,97]), fHeight=10)
# statTest(list(PCA.LDA.BV.I[, 97],
#               PCA.LDA.CV.I[,97],
#               PCA.LDA.BV.O[,97],
#               PCA.LDA.CV.O[,97],
#               PLS.LDA.BV.I[,97],
#               PLS.LDA.CV.I[,97],
#               PLS.LDA.BV.O[,97],
#               PLS.LDA.CV.O[,97]), item1=c(1,2,5,6), item2=c(3,4,7,8), colLine=rgb(0.5,0.5,1),colStar=rgb(0,0,1), fHeight=6)
# par(xaxt="n",yaxt="n",bty="n",mar=rep(0,4))
# plot(0:10/10,0:10/10,type="n",xlab="",ylab="")
legend('bottomright',legend=c("Test Result", "Val Result")
       ,lty=c(1,1,2,2),col=c(1,rgb(0.4,0.4,0.4)), pch=c(0,0),cex=0.8, bty='n')
dev.off()

PCA.SVML.BV.I <- c()
PCA.SVML.CV.I <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pca(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=svm, Pred.fn=predSVM, kernel='linear')
  PCA.SVML.BV.I <- rbind(PCA.SVML.BV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))

  TMP <- Make.classify.pca(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=svm, Pred.fn=predSVM, kernel='linear')
  PCA.SVML.CV.I <- rbind(PCA.SVML.CV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
}

# PLS+BV/CV Inside
PLS.SVML.BV.I <- c()
PLS.SVML.CV.I <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pls(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=svm, Pred.fn=predSVM, kernel='linear')
  PLS.SVML.BV.I <- rbind(PLS.SVML.BV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))

  TMP <- Make.classify.pls(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=svm, Pred.fn=predSVM, kernel='linear')
  PLS.SVML.CV.I <- rbind(PLS.SVML.CV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
}

# PCA+BV/CV Outside
PCA.SVML.BV.O <- c()
PCA.SVML.CV.O <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pca1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=svm, Pred.fn=predSVM, kernel='linear')
  PCA.SVML.BV.O <- rbind(PCA.SVML.BV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))

  TMP <- Make.classify.pca1(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=svm, Pred.fn=predSVM, kernel='linear')
  PCA.SVML.CV.O <- rbind(PCA.SVML.CV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
}

# PLS+BV/CV Outside
PLS.SVML.BV.O <- c()
PLS.SVML.CV.O <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pls1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=svm, Pred.fn=predSVM, kernel='linear')
  PLS.SVML.BV.O <- rbind(PLS.SVML.BV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))

  TMP <- Make.classify.pls1(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=svm, Pred.fn=predSVM, kernel='linear')
  PLS.SVML.CV.O <- rbind(PLS.SVML.CV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
}

save(file=paste(getwd(), "/Results/SVM_L_Mice.RData", sep=''), list=c("PCA.SVML.BV.I",
                                                                        "PCA.SVML.CV.I",
                                                                        "PLS.SVML.BV.I",
                                                                        "PLS.SVML.CV.I",
                                                                        "PCA.SVML.BV.O",
                                                                        "PCA.SVML.CV.O",
                                                                        "PLS.SVML.BV.O",
                                                                        "PLS.SVML.CV.O"))
load(file=paste(getwd(), "/Results/SVM_L_Mice.RData", sep=''))
k=1
pdf('Mice_SVM_L.pdf', width=7, height=3)
plot.accSd(Index=3:50, DATA=list(PCA.SVML.BV.I[k, 1:48],
                                 PCA.SVML.CV.I[k, 1:48],
                                 PLS.SVML.BV.I[k, 1:48],
                                 PLS.SVML.CV.I[k, 1:48],
                                 PCA.SVML.BV.O[k, 1:48],
                                 PCA.SVML.CV.O[k, 1:48],
                                 PLS.SVML.BV.O[k, 1:48],
                                 PLS.SVML.CV.O[k, 1:48]), 
           SD=list(PCA.SVML.BV.I[k, 49:96],
                   PCA.SVML.CV.I[k, 49:96],
                   PLS.SVML.BV.I[k, 49:96],
                   PLS.SVML.CV.I[k, 49:96],
                   PCA.SVML.BV.O[k, 49:96],
                   PCA.SVML.CV.O[k, 49:96],
                   PLS.SVML.BV.O[k, 49:96],
                   PLS.SVML.CV.O[k, 49:96]), 
           Names,YLIM=NULL,TYPE="l",MAIN="Mice Data (SVM.L)")
dev.off()

pdf('box_Mice_SVM_L.pdf', width=5, height=3)
# nf <- layout(matrix(c(1,2),nrow=1,byrow=FALSE), widths=c(4,1))
plot.Box(DATA=list(PCA.SVML.BV.I[, 97],
                   PCA.SVML.CV.I[,97],
                   PCA.SVML.BV.O[,97],
                   PCA.SVML.CV.O[,97],
                   PLS.SVML.BV.I[,97],
                   PLS.SVML.CV.I[,97],
                   PLS.SVML.BV.O[,97],
                   PLS.SVML.CV.O[,97]), 
         MAIN="Mice Data (SVM.L)", xtext=TRUE, method="", mar=c(3,4,2,2), cex.main=0.8)

plot.Box(DATA=list(apply(PCA.SVML.BV.I[,1:48],1,max),
                   apply(PCA.SVML.CV.I[,1:48],1,max),
                   apply(PCA.SVML.BV.O[,1:48],1,max),
                   apply(PCA.SVML.CV.O[,1:48],1,max),
                   apply(PLS.SVML.BV.I[,1:48],1,max),
                   apply(PLS.SVML.CV.I[,1:48],1,max),
                   apply(PLS.SVML.BV.O[,1:48],1,max),
                   apply(PLS.SVML.CV.O[,1:48],1,max)),
         MAIN="", xtext=FALSE, method="", mar=c(3,4,2,2), at=(1:8)+0.3, add=TRUE,border=rgb(0.4,0.4,0.4))
# statTest(list(PCA.SVML.BV.I[, 97],
#               PCA.SVML.CV.I[,97],
#               PCA.SVML.BV.O[,97],
#               PCA.SVML.CV.O[,97],
#               PLS.SVML.BV.I[,97],
#               PLS.SVML.CV.I[,97],
#               PLS.SVML.BV.O[,97],
#               PLS.SVML.CV.O[,97]), fHeight=10)
# statTest(list(PCA.SVML.BV.I[, 97],
#               PCA.SVML.CV.I[,97],
#               PCA.SVML.BV.O[,97],
#               PCA.SVML.CV.O[,97],
#               PLS.SVML.BV.I[,97],
#               PLS.SVML.CV.I[,97],
#               PLS.SVML.BV.O[,97],
#               PLS.SVML.CV.O[,97]), item1=c(1,2,5,6), item2=c(3,4,7,8), colLine=rgb(0.5,0.5,1),colStar=rgb(0,0,1), fHeight=6)
# par(xaxt="n",yaxt="n",bty="n",mar=rep(0,4))
# plot(0:10/10,0:10/10,type="n",xlab="",ylab="")
legend('bottomright',legend=c("Test Result", "Val Result")
       ,lty=c(1,1,2,2),col=c(1,rgb(0.4,0.4,0.4)), pch=c(0,0),cex=0.8, bty='n')
dev.off()

PCA.SVMR.BV.I <- c()
PCA.SVMR.CV.I <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pca(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=svm, Pred.fn=predSVM, kernel='radial')
  PCA.SVMR.BV.I <- rbind(PCA.SVMR.BV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))

  TMP <- Make.classify.pca(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=svm, Pred.fn=predSVM, kernel='radial')
  PCA.SVMR.CV.I <- rbind(PCA.SVMR.CV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
}

# PLS+BV/CV Inside
PLS.SVMR.BV.I <- c()
PLS.SVMR.CV.I <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pls(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=svm, Pred.fn=predSVM, kernel='radial')
  PLS.SVMR.BV.I <- rbind(PLS.SVMR.BV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))

  TMP <- Make.classify.pls(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                           Test=spec[ixTest,], testLabels=labels[ixTest],
                           Class.fn=svm, Pred.fn=predSVM, kernel='radial')
  PLS.SVMR.CV.I <- rbind(PLS.SVMR.CV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
}

# PCA+BV/CV Outside
PCA.SVMR.BV.O <- c()
PCA.SVMR.CV.O <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pca1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=svm, Pred.fn=predSVM, kernel='radial')
  PCA.SVMR.BV.O <- rbind(PCA.SVMR.BV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))

  TMP <- Make.classify.pca1(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=svm, Pred.fn=predSVM, kernel='radial')
  PCA.SVMR.CV.O <- rbind(PCA.SVMR.CV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
}

# PLS+BV/CV Outside
PLS.SVMR.BV.O <- c()
PLS.SVMR.CV.O <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pls1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=svm, Pred.fn=predSVM, kernel='radial')
  PLS.SVMR.BV.O <- rbind(PLS.SVMR.BV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))

  TMP <- Make.classify.pls1(Spec=spec[-ixTest, ], Labels=labels[-ixTest],
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=svm, Pred.fn=predSVM, kernel='radial')
  PLS.SVMR.CV.O <- rbind(PLS.SVMR.CV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
}

save(file=paste(getwd(), "/Results/SVM_R_Mice.RData", sep=''), list=c("PCA.SVMR.BV.I",
                                                                        "PCA.SVMR.CV.I",
                                                                        "PLS.SVMR.BV.I",
                                                                        "PLS.SVMR.CV.I",
                                                                        "PCA.SVMR.BV.O",
                                                                        "PCA.SVMR.CV.O",
                                                                        "PLS.SVMR.BV.O",
                                                                        "PLS.SVMR.CV.O"))

load(file=paste(getwd(), "/Results/SVM_R_Mice.RData", sep=''))
k=1
pdf('Mice_SVM_R.pdf', width=7, height=3)
plot.accSd(Index=3:50, DATA=list(PCA.SVMR.BV.I[k, 1:48],
                                 PCA.SVMR.CV.I[k, 1:48],
                                 PLS.SVMR.BV.I[k, 1:48],
                                 PLS.SVMR.CV.I[k, 1:48],
                                 PCA.SVMR.BV.O[k, 1:48],
                                 PCA.SVMR.CV.O[k, 1:48],
                                 PLS.SVMR.BV.O[k, 1:48],
                                 PLS.SVMR.CV.O[k, 1:48]), 
           SD=list(PCA.SVMR.BV.I[k, 49:96],
                   PCA.SVMR.CV.I[k, 49:96],
                   PLS.SVMR.BV.I[k, 49:96],
                   PLS.SVMR.CV.I[k, 49:96],
                   PCA.SVMR.BV.O[k, 49:96],
                   PCA.SVMR.CV.O[k, 49:96],
                   PLS.SVMR.BV.O[k, 49:96],
                   PLS.SVMR.CV.O[k, 49:96]), 
           Names,YLIM=NULL,TYPE="l",MAIN="Mice Data (SVM.R)")
dev.off()

pdf('box_Mice_SVM_R.pdf', width=5, height=3)
# nf <- layout(matrix(c(1,2),nrow=1,byrow=FALSE), widths=c(4,1))
plot.Box(DATA=list(PCA.SVMR.BV.I[, 97],
                   PCA.SVMR.CV.I[,97],
                   PCA.SVMR.BV.O[,97],
                   PCA.SVMR.CV.O[,97],
                   PLS.SVMR.BV.I[,97],
                   PLS.SVMR.CV.I[,97],
                   PLS.SVMR.BV.O[,97],
                   PLS.SVMR.CV.O[,97]), 
         MAIN="Mice Data (SVM.R)", xtext=TRUE, method="", mar=c(3,4,2,2), cex.main=0.8)

plot.Box(DATA=list(apply(PCA.SVMR.BV.I[,1:48],1,max),
                   apply(PCA.SVMR.CV.I[,1:48],1,max),
                   apply(PCA.SVMR.BV.O[,1:48],1,max),
                   apply(PCA.SVMR.CV.O[,1:48],1,max),
                   apply(PLS.SVMR.BV.I[,1:48],1,max),
                   apply(PLS.SVMR.CV.I[,1:48],1,max),
                   apply(PLS.SVMR.BV.O[,1:48],1,max),
                   apply(PLS.SVMR.CV.O[,1:48],1,max)),
         MAIN="", xtext=FALSE, method="", mar=c(3,4,2,2), at=(1:8)+0.3, add=TRUE,border=rgb(0.4,0.4,0.4))
# statTest(list(PCA.SVMR.BV.I[, 97],
#               PCA.SVMR.CV.I[,97],
#               PCA.SVMR.BV.O[,97],
#               PCA.SVMR.CV.O[,97],
#               PLS.SVMR.BV.I[,97],
#               PLS.SVMR.CV.I[,97],
#               PLS.SVMR.BV.O[,97],
#               PLS.SVMR.CV.O[,97]), fHeight=10)
# statTest(list(PCA.SVMR.BV.I[, 97],
#               PCA.SVMR.CV.I[,97],
#               PCA.SVMR.BV.O[,97],
#               PCA.SVMR.CV.O[,97],
#               PLS.SVMR.BV.I[,97],
#               PLS.SVMR.CV.I[,97],
#               PLS.SVMR.BV.O[,97],
#               PLS.SVMR.CV.O[,97]), item1=c(1,2,5,6), item2=c(3,4,7,8), colLine=rgb(0.5,0.5,1),colStar=rgb(0,0,1), fHeight=6)
# par(xaxt="n",yaxt="n",bty="n",mar=rep(0,4))
# plot(0:10/10,0:10/10,type="n",xlab="",ylab="")
legend('bottomright',legend=c("Test Result", "Val Result")
       ,lty=c(1,1,2,2),col=c(1,rgb(0.4,0.4,0.4)), pch=c(0,0),cex=0.8, bty='n')
dev.off()

PCA.RF.BV.I <- c()
PCA.RF.CV.I <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pca(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest], 
                           Test=spec[ixTest,], testLabels=labels[ixTest], 
                           Class.fn=randomForest, Pred.fn=predRF)
  PCA.RF.BV.I <- rbind(PCA.RF.BV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
  
  TMP <- Make.classify.pca(Spec=spec[-ixTest, ], Labels=labels[-ixTest], 
                           Test=spec[ixTest,], testLabels=labels[ixTest], 
                           Class.fn=randomForest, Pred.fn=predRF)
  PCA.RF.CV.I <- rbind(PCA.RF.CV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
}

# PLS+BV/CV Inside
PLS.RF.BV.I <- c()
PLS.RF.CV.I <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pls(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest], 
                           Test=spec[ixTest,], testLabels=labels[ixTest], 
                           Class.fn=randomForest, Pred.fn=predRF)
  PLS.RF.BV.I <- rbind(PLS.RF.BV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
  
  TMP <- Make.classify.pls(Spec=spec[-ixTest, ], Labels=labels[-ixTest], 
                           Test=spec[ixTest,], testLabels=labels[ixTest], 
                           Class.fn=randomForest, Pred.fn=predRF)
  PLS.RF.CV.I <- rbind(PLS.RF.CV.I, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
}

# PCA+BV/CV Outside
PCA.RF.BV.O <- c()
PCA.RF.CV.O <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pca1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest], 
                            Test=spec[ixTest,], testLabels=labels[ixTest],
                            Class.fn=randomForest, Pred.fn=predRF)
  PCA.RF.BV.O <- rbind(PCA.RF.BV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
  
  TMP <- Make.classify.pca1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], 
                            Test=spec[ixTest,], testLabels=labels[ixTest], 
                            Class.fn=randomForest, Pred.fn=predRF)
  PCA.RF.CV.O <- rbind(PCA.RF.CV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optPC))
}

# PLS+BV/CV Outside
PLS.RF.BV.O <- c()
PLS.RF.CV.O <- c()
for(i in 1:length(fCV))
{
  ixTest <- fCV[[i]]
  TMP <- Make.classify.pls1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], Batch=batches[-ixTest], 
                            Test=spec[ixTest,], testLabels=labels[ixTest], 
                            Class.fn=randomForest, Pred.fn=predRF)
  PLS.RF.BV.O <- rbind(PLS.RF.BV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
  
  TMP <- Make.classify.pls1(Spec=spec[-ixTest, ], Labels=labels[-ixTest], 
                            Test=spec[ixTest,], testLabels=labels[ixTest], 
                            Class.fn=randomForest, Pred.fn=predRF)
  PLS.RF.CV.O <- rbind(PLS.RF.CV.O, c(rowMeans(TMP$SENS), apply(TMP$SENS,1,sd), TMP$testSens, TMP$optCp))
}

save(file=paste(getwd(), "/Results/RF_Mice.RData", sep=''), list=c("PCA.RF.BV.I",
                                                                       "PCA.RF.CV.I",
                                                                       "PLS.RF.BV.I",
                                                                       "PLS.RF.CV.I",
                                                                       "PCA.RF.BV.O",
                                                                       "PCA.RF.CV.O",
                                                                       "PLS.RF.BV.O",
                                                                       "PLS.RF.CV.O"))

load(file=paste(getwd(), "/Results/RF_Mice.RData", sep=''))
k=1
pdf('Mice_RanFor.pdf', width=7, height=3)
plot.accSd(Index=3:50, DATA=list(PCA.RF.BV.I[k, 1:48],
                                 PCA.RF.CV.I[k, 1:48],
                                 PLS.RF.BV.I[k, 1:48],
                                 PLS.RF.CV.I[k, 1:48],
                                 PCA.RF.BV.O[k, 1:48],
                                 PCA.RF.CV.O[k, 1:48],
                                 PLS.RF.BV.O[k, 1:48],
                                 PLS.RF.CV.O[k, 1:48]), 
           SD=list(PCA.RF.BV.I[k, 49:96],
                   PCA.RF.CV.I[k, 49:96],
                   PLS.RF.BV.I[k, 49:96],
                   PLS.RF.CV.I[k, 49:96],
                   PCA.RF.BV.O[k, 49:96],
                   PCA.RF.CV.O[k, 49:96],
                   PLS.RF.BV.O[k, 49:96],
                   PLS.RF.CV.O[k, 49:96]), 
           Names,YLIM=NULL,TYPE="l",MAIN="Mice Data (RanFor)")
dev.off()

pdf('box_Mice_RanFor.pdf', width=5, height=3)
# nf <- layout(matrix(c(1,2),nrow=1,byrow=FALSE), widths=c(4,1))
plot.Box(DATA=list(PCA.RF.BV.I[, 97],
                   PCA.RF.CV.I[,97],
                   PCA.RF.BV.O[,97],
                   PCA.RF.CV.O[,97],
                   PLS.RF.BV.I[,97],
                   PLS.RF.CV.I[,97],
                   PLS.RF.BV.O[,97],
                   PLS.RF.CV.O[,97]), 
         MAIN="Mice Data (RanFor)", xtext=TRUE, method="", mar=c(3,4,2,2), cex.main=0.8)

plot.Box(DATA=list(apply(PCA.RF.BV.I[,1:48],1,max),
                   apply(PCA.RF.CV.I[,1:48],1,max),
                   apply(PCA.RF.BV.O[,1:48],1,max),
                   apply(PCA.RF.CV.O[,1:48],1,max),
                   apply(PLS.RF.BV.I[,1:48],1,max),
                   apply(PLS.RF.CV.I[,1:48],1,max),
                   apply(PLS.RF.BV.O[,1:48],1,max),
                   apply(PLS.RF.CV.O[,1:48],1,max)),
         MAIN="", xtext=FALSE, method="", mar=c(3,4,2,2), at=(1:8)+0.3, add=TRUE,border=rgb(0.4,0.4,0.4))
# statTest(list(PCA.RF.BV.I[, 97],
#               PCA.RF.CV.I[,97],
#               PCA.RF.BV.O[,97],
#               PCA.RF.CV.O[,97],
#               PLS.RF.BV.I[,97],
#               PLS.RF.CV.I[,97],
#               PLS.RF.BV.O[,97],
#               PLS.RF.CV.O[,97]), fHeight=10)
# statTest(list(PCA.RF.BV.I[, 97],
#               PCA.RF.CV.I[,97],
#               PCA.RF.BV.O[,97],
#               PCA.RF.CV.O[,97],
#               PLS.RF.BV.I[,97],
#               PLS.RF.CV.I[,97],
#               PLS.RF.BV.O[,97],
#               PLS.RF.CV.O[,97]), item1=c(1,2,5,6), item2=c(3,4,7,8), colLine=rgb(0.5,0.5,1),colStar=rgb(0,0,1), fHeight=6)
# par(xaxt="n",yaxt="n",bty="n",mar=rep(0,4))
# plot(0:10/10,0:10/10,type="n",xlab="",ylab="")
legend('bottomright',legend=c("Test Result", "Val Result")
       ,lty=c(1,1,2,2),col=c(1,rgb(0.4,0.4,0.4)), pch=c(0,0),cex=0.8, bty='n')
dev.off()