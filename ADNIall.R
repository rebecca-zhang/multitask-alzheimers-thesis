# Total ADNI 3T baseline data
# Rebecca Zhang, 2015

setwd("~/Dropbox/Multitask Thesis")
data1 = read.csv('UCSFFSX51_ADNI1_3T_08_01_14.csv',header=TRUE)
data2 = read.csv('UCSFFSX51_08_01_14.csv',header=TRUE)

subjects1 = sort(unique(data1$RID[data1$VISCODE == "bl"]))
subjects2 = sort(unique(data2$RID[data2$VISCODE2 == "scmri"]))
# check for overlap
sum(subjects1 %in% subjects2)
allsubjects = sort(c(subjects1,subjects2))
data1 = data1[which(data1$RID %in% allsubjects),]
data2 = data2[which(data2$RID %in% allsubjects),]
data1.info = data1[,which(names(data1) %in% c("RID","VISCODE","OVERALLQC"))]
data1 = cbind(data1.info,data1[,-(1:20)])
data2.info = data2[,which(names(data2) %in% c("RID","VISCODE2","OVERALLQC"))]
names(data2.info) = c("RID","VISCODE","OVERALLQC")
data2 = cbind(data2.info,data2[,-(1:23)])
deletedFeatures = c('ST100SV','ST122SV','ST126SV','ST22CV','ST22SA','ST22TA','ST22TS',
                    'ST28CV','ST33SV','ST41SV','ST63SV','ST67SV','ST81CV','ST81SA','ST81TA',
                    'ST81TS','ST87CV','ST92SV','ST8SV') 
deletedFeatures %in% colnames(data1)
deletedFeatures %in% colnames(data2)
data1 = data1[,-which(colnames(data1) %in% deletedFeatures)]
data2 = data2[,-which(colnames(data2) %in% deletedFeatures)]
# use data2 features, fewer
data1 = data1[,which(colnames(data1) %in% colnames(data2))]
data1 = data1[,match(names(data1),names(data2))]
data = rbind(data1,data2)
# baseline
baseline = which(data$VISCODE == 'bl')
baseline = c(baseline, which(data$VISCODE == 'scmri'))
data = data[baseline,]
# quality control
data = data[-which(data$OVERALLQC == "Fail"),]
data = data[-which(data$OVERALLQC == "Hippocampus Only"),]
subjects.order = unique(data$RID)
subjects.n = length(subjects.order)
# check missing features
max(colSums(is.na(data)))
features = data[,-c(1,2,3,344)]
features.order = names(features)
features.n = length(features.order)

# read in target vals
MMSEdata = read.csv('MMSE.csv',header=TRUE)
MMSEdata = MMSEdata[which(MMSEdata$RID %in% subjects.order),]
ADASdata = read.csv('ADAS_ADNIGO2.csv',header=TRUE)
ADASdata = ADASdata[which(ADASdata$RID %in% subjects.order),]
MMSEtimes = c("sc","m06","m12","m24","m36","m48")
ADAStimes = c("bl","m06","m12","m24","m36","m48")
times.n = length(MMSEtimes)
# remove anything not in the wanted times
MMSEdata = MMSEdata[which(MMSEdata$VISCODE2 %in% MMSEtimes),]
ADASdata = ADASdata[which(ADASdata$VISCODE2 %in% ADAStimes),]
# remove NA or negative records
MMSEdata = MMSEdata[-which(is.na(MMSEdata$MMSCORE)),]
MMSEdata = MMSEdata[-which(MMSEdata$MMSCORE < 0),]
ADASdata = ADASdata[-which(is.na(ADASdata$TOTAL13)),]
ADASdata = ADASdata[which(ADASdata$TOTAL13 >= 0),]
length(unique(MMSEdata$RID)) # check these are equal to subjects.n, or else have prob
length(unique(ADASdata$RID))

# STORE DATA
X = matrix(NA,subjects.n,features.n)
for (i in 1:dim(data)[1])
{
  #print(i)
  subjects.index = match(data$RID[i],subjects.order)
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

# read in MMSE
MMSE = matrix(NA,subjects.n,times.n)
for (i in 1:dim(MMSEdata)[1]) 
{
  times.index = match(MMSEdata$VISCODE2[i],MMSEtimes)
  subjects.index = match(MMSEdata$RID[i],subjects.order)
  MMSE[subjects.index,times.index] = MMSEdata$MMSCORE[i]
  print(i)
}
colnames(MMSE) = MMSEtimes
colSums(is.na(MMSE)==FALSE)

# read in ADAS
ADAS = matrix(NA,subjects.n,times.n)
for (i in 1:dim(ADASdata)[1])
{
  times.index = match(ADASdata$VISCODE2[i],ADAStimes)
  subjects.index = match(ADASdata$RID[i],subjects.order)
  ADAS[subjects.index,times.index] = ADASdata$TOTAL13[i]
  print(i)
}
colnames(ADAS) = ADAStimes
colSums(is.na(ADAS)==FALSE)

# add baseline MMSE to data matrix
X = cbind(X,MMSE[,1])
# remove baseline
MMSE = MMSE[,-1]
ADAS = ADAS[,-1]

write.csv(X,file="ADNI2baseline.csv",row.names=F,na="-4")
write.csv(MMSE,file="ADNI2MMSE.csv",row.names=F,na="-4")
write.csv(ADAS,file="ADNI2ADAS.csv",row.names=F,na="-4")

#-----------------------TRAINING AND PREDICTION-------------------------------
set.seed(123)
#param = exp(seq(log(3),log(5000),0.2))
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
trainset = matrix(0,floor(subjects.n*0.9),20)
for (cut in 1:20)
{
  trainset[,cut] = sort(sample(1:subjects.n,floor(subjects.n * 0.9)))
}
write.table(trainset,file='ADNI2trainindices.csv',sep=",",col.names=F,row.names=F)
for (cut in 1:20)
{
  train.indices = trainset[,cut]
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


# -------------------------------------------MODELS--------------------------------------------
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
    #print(i)
    Yi = Y[which(is.na(Y[,i])==FALSE),i]
    Xi = X[which(is.na(Y[,i])==FALSE),]
    cv.model = cv.glmnet(x=Xi,y=Yi,family="gaussian",alpha=1,nfolds=5,standardize=F)
    lambda = cv.model$lambda.min
    W[,i] = coef(cv.model,s="lambda.min")[-1] #don't use intercept, should be 0 anyway
  }
  return(W)
}


# ------------------------------------Performance Indicators -------------------------
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

