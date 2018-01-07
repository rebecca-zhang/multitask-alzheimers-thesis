# ADNI 1 replication
# Rebecca Zhang, 2015

# All input is 1.5 T data

library(zoo)
library(abind)
library(glmnet)

setwd("~/Dropbox/Multitask Thesis")
MRIdata = read.csv('UCSFFSX_08_01_14.csv',header=TRUE)
# check the right features are removed
deletedFeatures = c('ST100SV','ST122SV','ST126SV','ST22CV','ST22SA','ST22TA','ST22TS',
                    'ST28CV','ST33SV','ST41SV','ST63SV','ST67SV','ST81CV','ST81SA','ST81TA',
                    'ST81TS','ST87CV','ST92SV','ST8SV') 
deletedFeatures %in% colnames(MRIdata)
MRIdata = MRIdata[,-which(colnames(MRIdata) %in% deletedFeatures)]
# remove irrelevant columns
MRIdata.info = MRIdata[,which(names(MRIdata) %in% c("RID","VISCODE","OVERALLQC"))]
MRIdata = cbind(MRIdata.info,MRIdata[,-(1:19)])

# eliminate image records with failed overall QC
#MRIdata = MRIdata[which(MRIdata$OVERALLQC == "Pass"),]
MRIdata = MRIdata[-which(MRIdata$OVERALLQC == "Fail"),]
MRIdata = MRIdata[-which(MRIdata$OVERALLQC == "Hippocampus Only"),]

# include only baseline
baseline = MRIdata[which(MRIdata$VISCODE=='sc'),]

subjects.order = sort(unique(baseline$RID))
subjects.n = length(subjects.order)

# Y target scores
MMSEdata = read.csv('MMSE.csv',header=TRUE)
MMSEdata = MMSEdata[which(MMSEdata$Phase == "ADNI1"),]
MMSEdata = MMSEdata[which(MMSEdata$RID %in% subjects.order),]
ADASdata = read.csv('ADASSCORES.csv',header=TRUE)
ADASdata = ADASdata[which(ADASdata$RID %in% subjects.order),]
MMSEtimes = c("sc","m06","m12","m24","m36","m48")
ADAStimes = c("bl","m06","m12","m24","m36","m48")
times.n = length(MMSEtimes)
# remove anything not in the wanted times
MMSEdata = MMSEdata[which(MMSEdata$VISCODE2 %in% MMSEtimes),]
ADASdata = ADASdata[which(ADASdata$VISCODE %in% ADAStimes),]
# remove those without score
MMSEdata = MMSEdata[-which(MMSEdata$MMSCORE < 0),]
ADASdata = ADASdata[which(ADASdata$TOTAL11 >= 0),]
# update subjects to those who have m06 score
length(unique(MMSEdata$RID[which(MMSEdata$VISCODE2=="m06")]))
length(unique(ADASdata$RID[which(ADASdata$VISCODE=="m06")])) # make sure these equal
subjects.order = sort(unique(MMSEdata$RID[which(MMSEdata$VISCODE2=="m06")]))
subjects.n = length(subjects.order)
# update MRI, scores accordingly
baseline = baseline[which(baseline$RID %in% subjects.order),]
features = baseline[,-(1:3)]
features = features[,-dim(features)[2]] # get rid of update_stamp
features.order = colnames(features)
features.n = length(features.order)
MMSEdata = MMSEdata[which(MMSEdata$RID %in% subjects.order),]
ADASdata = ADASdata[which(ADASdata$RID %in% subjects.order),]

# STORE DATA
X = matrix(NA,subjects.n,features.n)
for (i in 1:dim(baseline)[1])
{
  #print(i)
  subjects.index = match(baseline$RID[i],subjects.order)
  for(j in 1:features.n)
    X[subjects.index,j] = features[i,j]
}
max(colSums(is.na(X)))
missing = which(colSums(is.na(X))>0)
# replace with avg value
for(i in 1:length(missing))
{
  feature = missing[i]
  X[is.na(X[,feature]),feature] = mean(X[,feature],na.rm=T)
}
max(colSums(is.na(X)))

# read in scores
MMSE = matrix(NA,subjects.n,times.n)
for (i in 1:dim(MMSEdata)[1]) 
{
  times.index = match(MMSEdata$VISCODE2[i],MMSEtimes)
  subjects.index = match(MMSEdata$RID[i],subjects.order)
  MMSE[subjects.index,times.index] = MMSEdata$MMSCORE[i]
  #print(i)
}
colnames(MMSE) = MMSEtimes

# read in ADAS data
ADAS = matrix(NA, subjects.n, times.n)
for (i in 1:dim(ADASdata)[1])
{
  times.index = match(ADASdata$VISCODE[i],ADAStimes)
  subjects.index = match(ADASdata$RID[i],subjects.order)
  ADAS[subjects.index,times.index] = ADASdata$TOTAL11[i]
}
colnames(ADAS) = ADAStimes
colSums(is.na(MMSE)==FALSE)
colSums(is.na(ADAS)==FALSE)


# add baseline MMSE to data matrix
X = cbind(X,MMSE[,1])
# remove baseline
MMSE = MMSE[,-1]
# remove bl
ADAS = ADAS[,-1]

write.csv(X,file="ADNI1baseline.csv",row.names=F,na="-4")
write.csv(MMSE,file="ADNI1MMSE.csv",row.names=F,na="-4")
write.csv(ADAS,file="ADNI1ADAS.csv",row.names=F,na="-4")

# ------------- TRAINING AND PREDICTION -------------------------------------
set.seed(123)
#param = exp(seq(log(3),log(5000),0.2))
trainindices = read.csv('ADNI1trainindices.csv',header=F)
nMSE.ridge.MMSE = c()
nMSE.lasso.MMSE = c()
wR.ridge.MMSE = c()
wR.lasso.MMSE = c()
rMSE.ridge.MMSE = c()
rMSE.lasso.MMSE = c()
nMSE.ridge.ADAS = c()
nMSE.lasso.ADAS = c()
wR.ridge.ADAS = c()
wR.lasso.ADAS = c()
rMSE.ridge.ADAS = c()
rMSE.lasso.ADAS = c()
for (cut in 1:20)
{
  #train.indices = sort(sample(1:subjects.n,floor(subjects.n * 0.9)))
  train.indices = trainindices[,cut]
  X.train = X[train.indices,] 
  MMSE.train = MMSE[train.indices,]
  ADAS.train = ADAS[train.indices,]
  Xmean = colMeans(X.train)
  Xsd = c()
  for (j in 1:dim(X.train)[2])
  {
    Xsd = c(Xsd,sd(X.train[,j]))
    X.train[,j] = (X.train[,j] - mean(X.train[,j]))/sd(X.train[,j])
  }
  #XMMSE.train = sweep(sweep(XMMSE.train,2,colMeans(XMMSE.train)),2,apply(XMMSE.train,2,sd),FUN="/")
  MMSEmean = colMeans(MMSE.train,na.rm=T)
  MMSEsd = c()
  for (j in 1:dim(MMSE.train)[2])
  {
    MMSEsd = c(MMSEsd,sd(MMSE.train[,j],na.rm=T))
    MMSE.train[,j] = (MMSE.train[,j] - mean(MMSE.train[,j],na.rm=T))/sd(MMSE.train[,j],na.rm=T)
  }
  ADASmean = colMeans(ADAS.train,na.rm=T)
  ADASsd = c()
  for (j in 1:dim(ADAS.train)[2])
  {
    ADASsd = c(ADASsd,sd(ADAS.train[,j],na.rm=T))
    ADAS.train[,j] = (ADAS.train[,j] - mean(ADAS.train[,j],na.rm=T))/sd(ADAS.train[,j],na.rm=T)
  }
  
  X.test = X[-train.indices,]
  MMSE.test = MMSE[-train.indices,]
  ADAS.test = ADAS[-train.indices,]
  for (j in 1:dim(X.test)[2])
    #X.test[,j] = (X.test[,j] - mean(X.test[,j]))/sd(X.test[,j])
    X.test[,j] = (X.test[,j] - Xmean[j])/Xsd[j]
  #XMMSE.test = sweep(sweep(XMMSE.test,2,colMeans(XMMSE.test)),2,apply(XMMSE.test,2,sd),FUN="/")
  #for (j in 1:dim(MMSE.test)[2])
  # MMSE.test[,j] = (MMSE.test[,j] - mean(MMSE.test[,j],na.rm=T))/sd(MMSE.test[,j],na.rm=T)
  
  predictedMMSEridge = X.test %*% ridge(X.train,MMSE.train)
  predictedMMSElasso = X.test %*% lasso(X.train,MMSE.train)
  predictedADASridge = X.test %*% ridge(X.train,ADAS.train)
  predictedADASlasso = X.test %*% lasso(X.train,ADAS.train)
  for (j in 1:dim(predictedMMSEridge)[2])
  {
    predictedMMSEridge[,j] = predictedMMSEridge[,j] * MMSEsd[j] + MMSEmean[j]
    predictedMMSElasso[,j] = predictedMMSElasso[,j] * MMSEsd[j] + MMSEmean[j]
  }
  for (j in 1:dim(predictedADASridge)[2])
  {
    predictedADASridge[,j] = predictedADASridge[,j] * ADASsd[j] + ADASmean[j]
    predictedADASlasso[,j] = predictedADASlasso[,j] * ADASsd[j] + ADASmean[j]
  }
  nMSE.ridge.MMSE = c(nMSE.ridge.MMSE, normalizedMeanSquared(MMSE.test,predictedMMSEridge))
  nMSE.lasso.MMSE = c(nMSE.lasso.MMSE, normalizedMeanSquared(MMSE.test,predictedMMSElasso))
  wR.ridge.MMSE = c(wR.ridge.MMSE, weightedR(MMSE.test,predictedMMSEridge))
  wR.lasso.MMSE = c(wR.lasso.MMSE, weightedR(MMSE.test,predictedMMSElasso))
  rMSE.ridge.MMSE = rbind(rMSE.ridge.MMSE, rootMeanSquared(MMSE.test,predictedMMSEridge))
  rMSE.lasso.MMSE = rbind(rMSE.lasso.MMSE, rootMeanSquared(MMSE.test,predictedMMSElasso))
  
  nMSE.ridge.ADAS = c(nMSE.ridge.ADAS, normalizedMeanSquared(ADAS.test,predictedADASridge))
  nMSE.lasso.ADAS = c(nMSE.lasso.ADAS, normalizedMeanSquared(ADAS.test,predictedADASlasso))
  wR.ridge.ADAS = c(wR.ridge.ADAS, weightedR(ADAS.test,predictedADASridge))
  wR.lasso.ADAS = c(wR.lasso.ADAS, weightedR(ADAS.test,predictedADASlasso))
  rMSE.ridge.ADAS = rbind(rMSE.ridge.ADAS, rootMeanSquared(ADAS.test,predictedADASridge))
  rMSE.lasso.ADAS = rbind(rMSE.lasso.ADAS, rootMeanSquared(ADAS.test,predictedADASlasso))
}
mean(nMSE.ridge.MMSE)
mean(wR.ridge.MMSE)
colMeans(rMSE.ridge.MMSE)
sd(nMSE.ridge.MMSE)
sd(wR.ridge.MMSE)
apply(rMSE.ridge.MMSE,2,sd)
mean(nMSE.ridge.ADAS)
mean(wR.ridge.ADAS)
colMeans(rMSE.ridge.ADAS)
sd(nMSE.ridge.ADAS)
sd(wR.ridge.ADAS)
apply(rMSE.ridge.ADAS,2,sd)

mean(nMSE.lasso.MMSE)
mean(wR.lasso.MMSE)
colMeans(rMSE.lasso.MMSE)
sd(nMSE.lasso.MMSE)
sd(wR.lasso.MMSE)
apply(rMSE.lasso.MMSE,2,sd)
mean(nMSE.lasso.ADAS)
mean(wR.lasso.ADAS)
colMeans(rMSE.lasso.ADAS)
sd(nMSE.lasso.ADAS)
sd(wR.lasso.ADAS)
apply(rMSE.lasso.ADAS,2,sd)

plot(MMSE.test,predictedMMSEridgeCV,xlim=c(10,30),ylim=c(10,30))
abline(0,1)
plot(MMSE.test,predictedMMSElassoCV,xlim=c(10,30),ylim=c(10,30))
abline(0,1)

minerror = min(avgspliterror)
bestparam = param[which(avgspliterror == min(avgspliterror))]




# ----------------------------------MODELS-----------------------------------------
# ridge, assumes X is standardized
ridge<-function(X,Y)
{
  d = dim(X)[2]
  n = colSums(is.na(Y)==FALSE)
  t = dim(Y)[2]
  W = matrix(NA,d,t)
  for (i in 1:t)
  {
    Yi = Y[which(is.na(Y[,i])==FALSE),i]
    Xi = X[which(is.na(Y[,i])==FALSE),]
    #W[,i] = solve(t(Xi)%*%Xi + theta*diag(d)) %*% t(Xi) %*% Yi
    #W[,i] = solve(t(Xi)%*%Xi + theta*diag(d),t(Xi) %*% Yi)
    cv.model = cv.glmnet(x=Xi,y=Yi,family="gaussian",alpha=0,nfolds=5,standardize=F)
    lambda = cv.model$lambda.min
    W[,i] = coef(cv.model,s="lambda.min")[-1]
  }
  return(W)
}

# lasso, assumes X is standardized
lasso<-function(X,Y)
{
  d = dim(X)[2]
  n = colSums(is.na(Y)==FALSE)
  t = dim(Y)[2]
  W = matrix(NA,d,t)
  for(i in 1:t)
  {
    print(i)
    Yi = Y[which(is.na(Y[,i])==FALSE),i]
    Xi = X[which(is.na(Y[,i])==FALSE),]
    cv.model = cv.glmnet(x=Xi,y=Yi,family="gaussian",alpha=1,nfolds=5,standardize=F)
    lambda = cv.model$lambda.min
    W[,i] = coef(cv.model,s="lambda.min")[-1] #don't use intercept, should be 0 anyway
  }
  return(W)
}

# ---------------------------------- Performance indicators --------------------------
rootMeanSquared <- function(Y, predicted)
{
  t = dim(Y)[2]
  rMSE = c()
  for(i in 1:t)
  {
    Yi = Y[which(is.na(Y[,i])==FALSE),i]
    Pi = predicted[which(is.na(Y[,i])==FALSE),i]
    ni = length(Yi)
    rMSE = c(rMSE,sqrt(sum((Yi-Pi)^2)/ni))
  }
  return(rMSE)
}
normalizedMeanSquared <- function(Y, predicted)
{
  t = dim(Y)[2]
  nMSE = 0
  n = 0
  for(i in 1:t)
  {
    Yi = Y[which(is.na(Y[,i])==FALSE),i]
    Pi = predicted[which(is.na(Y[,i])==FALSE),i]
    ni = length(Yi)
    if (ni > 1) # gives us NA otherwise
    {
      nMSE = nMSE + sum((Yi-Pi)^2)/var(Yi)
      n = n + ni
    }
  }
  nMSE = nMSE/n
  return(nMSE)
}
weightedR <- function(Y, predicted)
{
  t = dim(Y)[2]
  wR = 0
  n = 0
  for(i in 1:t)
  {
    Yi = Y[which(is.na(Y[,i])==FALSE),i]
    Pi = predicted[which(is.na(Y[,i])==FALSE),i]
    ni = length(Yi)
    if (ni > 2) # 2 or fewer points gives perfect correlation always
    {
      wR = wR + cor(Yi,Pi)*ni
      n = n + ni
    }
  } 
  wR = wR/n
  return(wR)
}
