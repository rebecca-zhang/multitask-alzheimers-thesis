# ADNI(1,GO,2) DATA INPUT
# Rebecca Zhang
# 3/15

library(zoo)
library(abind)

# X (feature matrix)
setwd("~/Dropbox/Multitask Thesis")
MRIdata = read.csv('UCSFFSX51_08_01_14.csv',header=TRUE)
MRIdata = MRIdata[MRIdata$IMAGETYPE == 'Non-Accelerated T1',] # remove accelerated T1 
MRIdata = MRIdata[MRIdata$VISCODE2 != 'nv',] # remove the two people who are NV
# check the right features are removed
c('ST100SV','ST122SV','ST126SV','ST22CV','ST22SA','ST22TA','ST22TS',
  'ST28CV','ST33SV','ST41SV','ST63SV','ST67SV','ST81CV','ST81SA','ST81TA',
  'ST81TS','ST87CV','ST92SV','ST8SV') %in% colnames(MRIdata)
MRIdata = MRIdata[,-which(colnames(MRIdata) == "ST8SV")]
# remove irrelevant columns
MRIdata.info = MRIdata[,which(names(MRIdata) %in% c("COLPROT","RID","VISCODE2","OVERALLQC"))]
colnames(MRIdata.info) = c("COLPROT","RID","VISCODE","OVERALLQC")
MRIdata = cbind(MRIdata.info,MRIdata[,-(1:23)])

# read in ADNI1 data
MRIdata1 = read.csv('UCSFFSX51_ADNI1_3T_08_01_14.csv',header=TRUE)
MRIdata1 = MRIdata1[,-which(is.na(MRIdata1[1,]))] # remove unused features
# replace all instances of bl with scmri
levels(MRIdata1$VISCODE)[match('bl',levels(MRIdata1$VISCODE))] = 'scmri'
MRIdata1$VISCODE[MRIdata1$VISCODE == 'bl'] = 'scmri'
# remove irrelevant columns
MRIdata1.info = MRIdata1[,which(names(MRIdata1) %in% c("RID","VISCODE","OVERALLQC"))]
MRIdata1.info = cbind(matrix("ADNI1",dim(MRIdata1)[1],1),MRIdata1.info)
colnames(MRIdata1.info) = c("COLPROT","RID","VISCODE","OVERALLQC")
MRIdata1 = cbind(MRIdata1.info,MRIdata1[,-(1:18)])
# check that features between two datasets are same
names(MRIdata1) == names(MRIdata)

# CONCATENATE TWO DATASETS
MRIdata = rbind(MRIdata,MRIdata1)
MRIdata$VISCODE = factor(MRIdata$VISCODE)

# eliminate patients w/o baseline MRI records
baselines = MRIdata$RID[MRIdata$VISCODE == 'scmri'] 
remove = which(MRIdata$RID %in% baselines == FALSE)
MRIdata = MRIdata[-remove,]

# eliminate image records with failed overall QC
MRIdata = MRIdata[-which(MRIdata$OVERALLQC == "Fail"),]

# eliminate subjects who only have 1 MRI reading
# onetime = names(which(table(MRIdata$RID) == 1))
# MRIdata = MRIdata[-which(MRIdata$RID %in% onetime),]

features = MRIdata[,-(1:4)]
features = features[,-dim(features)[2]] # get rid of update_stamp

times.order = sort(levels(MRIdata$VISCODE))
subjects.order = sort(unique(MRIdata$RID))
features.order = colnames(features)

times.n = length(times.order)
subjects.n = length(subjects.order)
features.n = length(features.order)
# move scmri, baseline, to the front
times.order = unlist(list(times.order[times.n],times.order[-times.n]))
X = array(NA,dim=c(subjects.n,features.n,times.n))
X.count = matrix(0,subjects.n,times.n)

# read in data
for (i in 1:dim(MRIdata)[1]) 
{
  times.index = which(times.order == MRIdata$VISCODE[i])
  subjects.index = which(subjects.order == MRIdata$RID[i])
  for (j in 1:features.n) # slower this way, looping, but can't figure out how to set at once??
  {
    if (X.count[subjects.index,times.index] > 0) # if done before
    {
      newavg = sum(X[subjects.index,j,times.index]*X.count[subjects.index,times.index], 
                   features[i,j],na.rm=TRUE)/(X.count[subjects.index,times.index] + 1)
      X[subjects.index,j,times.index] = newavg
    }
    else
    {
      X[subjects.index,j,times.index] = features[i,j]
    }
  }
  X.count[subjects.index,times.index] = X.count[subjects.index,times.index] + 1
}
# remove subjects with only 1 MRI
onlyone = which(rowSums(X.count==0) > 10)
X = X[-onlyone,,]
subjectsremoved = subjects.order[onlyone]
MRIdata = MRIdata[-which(MRIdata$RID %in% subjectsremoved),]
times.order = sort(levels(MRIdata$VISCODE))
subjects.order = sort(unique(MRIdata$RID))
features.order = colnames(features)
times.n = length(times.order)
subjects.n = length(subjects.order)
features.n = length(features.order)
# move scmri, baseline, to the front
times.order = unlist(list(times.order[times.n],times.order[-times.n]))


# Y (scores matrix)
Y = matrix(NA,subjects.n,times.n)
# scoreExists = matrix(0,subjects.n,times.n)
MMSEdata = read.csv('MMSE.csv',header=TRUE)
# include only subjects RIDs we also have MRI data for
MMSEdata = MMSEdata[which(MMSEdata$RID %in% subjects.order),]
# replace all instances of 'sc' with 'scmri'
levels(MMSEdata$VISCODE2)[match('sc',levels(MMSEdata$VISCODE2))] = 'scmri'
MMSEdata$VISCODE2[MMSEdata$VISCODE2 == 'sc'] = 'scmri'
# include only times for which we have MRIdata
MMSEdata = MMSEdata[MMSEdata$VISCODE2 %in% times.order,]
MMSEdata$VISCODE2 = factor(MMSEdata$VISCODE2, levels=times.order)
# sort subjects in order
MMSEdata = MMSEdata[order(MMSEdata$RID),]

# read in data
for (i in 1:dim(MMSEdata)[1]) 
{
  times.index = which(times.order == MMSEdata$VISCODE2[i])
  subjects.index = which(subjects.order == MMSEdata$RID[i])
  Y[subjects.index,times.index] = MMSEdata$MMSCORE[i]
  #print(i)
}
rownames(Y) = subjects.order
colnames(Y) = times.order

# complete missing features using avg value (?????)
# plot corresponding feature for all the time points for patients w/o missing data, 
# check if it varies over time, if not than the the avg
plot(1,1,xlim=c(1,12), ylim=range(X[,1,],na.rm=TRUE))
for (j in 1:(subjects.n-1000))
{
  lines(1:12, X.baseline[j,1,])
}

# interpolation over time
X.interpolated = X 
for (i in 1:dim(X.interpolated)[1])
{
  print(i)
  z = zoo(t(X.interpolated[i,,]))
  if(sum(colSums(is.na(z)) < 11) > 0) # more than just one MRI scan
  {
    if(sum(colSums(is.na(z)) < 11) < features.n) # not all features have 2+ entries
    {
      missingfeature = which(colSums(is.na(z)) >= 11)
      for (j in 1:features.n)
        if (j %in% missingfeature == FALSE)
          X.interpolated[i,j,] = na.approx(X.interpolated[i,j,],na.rm=FALSE)
    }
    else
      X.interpolated[i,,] = t(na.approx(z,na.rm=FALSE)) 
  }
}

Y.interpolated = Y
for (i in 1:dim(Y.interpolated)[1])
{
  print(i)
  y = zoo((Y.interpolated[i,]))
  if(sum(is.na(y)) <= 10)
    Y.interpolated[i,] = (na.approx(y,na.rm=FALSE))
}

# missing value counts
missingcount.X.interpolated = rep(0,times.n)
missingcount.Y.interpolated = rep(0,times.n)
for (i in 1:times.n)
{
  missingcount.X.interpolated[i] = sum(is.na(X.interpolated[,,i]))
  missingcount.Y.interpolated[i] = sum(is.na(Y.interpolated[,i]))
}
missingcount.X.interpolated
missingcount.Y.interpolated
# for baseline,m03,m06,m12 do the average to fill in missing values
for (i in 1:features.n)
{
  X.interpolated[is.na(X.interpolated[,i,1]),i,1] = mean(X.interpolated[,i,1],na.rm=TRUE)
  X.interpolated[is.na(X.interpolated[,i,2]),i,2] = mean(X.interpolated[,i,2],na.rm=TRUE)
  X.interpolated[is.na(X.interpolated[,i,3]),i,3] = mean(X.interpolated[,i,3],na.rm=TRUE)
  X.interpolated[is.na(X.interpolated[,i,4]),i,4] = mean(X.interpolated[,i,4],na.rm=TRUE)
}
# fill in MMSE scores for baseline,m03,m06,m12
Y.interpolated[is.na(Y.interpolated[,1]),1] = mean(Y.interpolated[,1],na.rm=TRUE)
Y.interpolated[is.na(Y.interpolated[,2]),2] = mean(Y.interpolated[,2],na.rm=TRUE)
Y.interpolated[is.na(Y.interpolated[,3]),3] = mean(Y.interpolated[,3],na.rm=TRUE)
Y.interpolated[is.na(Y.interpolated[,4]),4] = mean(Y.interpolated[,4],na.rm=TRUE)


# Eqn 3 - Analytical solution given in Appendix A
train3 <- function(theta1, theta2, X.3, Y.3)
{
  n = dim(X.3)[1]
  d = dim(X.3)[2]
  t = dim(Y.3)[2]
  H = matrix(0,t,t-1)
  for (i in 1:t)
    for (j in 1:t-1)
    {
      if(i==j) H[i,j] = 1
      if(i==j+1) H[i,j] = -1
    }
  # parameters theta need to be selected through CV
  sym1 = t(X.3)%*%X.3 + theta1*diag(d)
  sym2 = theta2*H%*%t(H)
  e1 = eigen(sym1)
  e2 = eigen(sym2)
  lambda1 = diag(e1$values)
  lambda2 = diag(e2$values)
  Q1 = e1$vectors
  Q2 = e2$vectors
  D = t(Q1) %*% t(X.3) %*% Y.3 %*% Q2
  What = D
  for (i in 1:dim(D)[1])
    for (j in 1:dim(D)[2])
      What[i,j] = What[i,j]/(e1$values[i] + e2$values[j])
  W = Q1 %*% What %*% t(Q2)
  return(W)
}

# Eqn 4 - Optimization problem. Fill in X, Y can be left missing.
#???????????? CHECK WITH HAN...the paper seems wrong
train4 <- function(theta1,theta2,X,Y) # X is 3-D
{
  d = dim(X)[2]
  t = dim(Y)[2]
  A = -theta2*diag(d)
  M = array(NA,dim=c(t,d,d))
  T = array(NA,dim=c(t,d,1))
  equationarray = matrix(0,t*d,t*d + 1)
  for(i in 1:t)
  {
    M[i,,] = theta1*diag(d) + 2*theta2*diag(d) + t(X[,,i])%*%X[,,i]
    T[i,,] = t(X[,,1])%*%Y[,i]
    if(i==1) {
      equationarray[1:d,1:(2*d)] = cbind(M[i,,],A)
    } else if(i==t) {
      equationarray[((t-1)*d+1):(t*d),((t-2)*d+1):(t*d)] = cbind(A,M[i,,])
    } else {
      equationarray[((i-1)*d+1):(i*d),((i-2)*d+1):((i+1)*d)] = cbind(A,M[i,,],A)
    }
  }
}

# Gradient descent
gradientDescent<-function(W0,alpha,theta1,theta2,S,X,Y,epsilon)
{
  W = W0
  H = matrix(0,t,t-1)
  for (i in 1:t)
    for (j in 1:t-1)
    {
      if(i==j) H[i,j] = 1
      if(i==j+1) H[i,j] = -1
    }
  D = 2*t(X)%*%(S * (X%*%W - Y))/dim(Y)[1] + 2*theta1*W + 2*theta2*W %*% H %*% t(H)
  i = 0
  while(norm(D,type="F") > epsilon)
  {
    #print(norm(D,type="F"))
    #print(norm(S* (X%*%W - Y),type="F")^2/dim(Y)[1] + theta1*norm(W,type="F")^2 + theta2*norm(W%*%H,type="F")^2)
    W = W - alpha*D
    D = 2*t(X)%*%(S * (X%*%W - Y))/dim(Y)[1] + 2*theta1*W + 2*theta2*W %*% H %*% t(H)
    i=i+1
  }
  return(W)
}

# train on the half and test on half (linear model prediction)
# For 3, X is the input data at the baseline that has not failed
train.indices = sort(sample(1:dim(X.interpolated)[1],floor(dim(X.interpolated)[1]/2)))
X.3.train = X.interpolated[train.indices,,1] #baseline only
Y.3.train = Y.interpolated[train.indices,1:4]
Y.3.train = sweep(sweep(Y.3.train,2,colMeans(Y.3.train)),2,apply(Y.3.train,2,sd),FUN="/")
X.3.train = sweep(sweep(X.3.train,2,colMeans(X.3.train)),2,apply(X.3.train,2,sd),FUN="/")

theta1.3 = seq(580, 600, 2)
theta2.3 = 400
rMSE.3 = matrix(NA,length(theta1.3),4)
nMSE.3 = rep(NA,length(theta1.3))
wR.3 = rep(0,length(theta1.3))
for (a in 1:length(theta1.3))
{
  coef.3 = train3(theta1.3[a],theta2.3,X.3.train,Y.3.train)
  X.3.test = X.interpolated[-train.indices,,1] 
  Y.3.test = Y.interpolated[-train.indices,1:4]
  Y.3.test = sweep(sweep(Y.3.test,2,colMeans(Y.3.test)),2,apply(Y.3.test,2,sd),FUN="/")
  X.3.test = sweep(sweep(X.3.test,2,colMeans(X.3.test)),2,apply(X.3.test,2,sd),FUN="/")
  predicted.3 = X.3.test %*% coef.3
  
  # calculate errors
  for(i in 1:4)
  {
    rMSE.3[a,i] = sqrt(sum((Y.3.test[,i]-predicted.3[,i])^2)/dim(Y.3.test)[1])
  }
  nMSE.3[a] = sum((Y.3.test - predicted.3)^2)/(dim(Y.3.test)[1]*dim(Y.3.test)[2])
  for (i in 1:4)
  {
    wR.3[a] = wR.3[a] + cor(Y.3.test[,i],predicted.3[,i])*dim(Y.3.test)[1]
  }
  wR.3[a] = wR.3[a] / (dim(Y.3.test)[1] * dim(Y.3.test)[2])
}

# 4
# X.4.train = X.interpolated[train.indices,,1:4] # check no missing values 
# Y.4.train = Y[train.indices,1:4]
# X.4.train = sweep(sweep(X.4.train,2,colMeans(X.4.train)),2,apply(X.4.train,2,sd),FUN="/")

# Gradient method
# with all target values
X.g.train = X.interpolated[train.indices,,1] #baseline only
Y.g.train = Y.interpolated[train.indices,1:4] # only these have all target vals
X.g.train = sweep(sweep(X.g.train,2,colMeans(X.g.train)),2,apply(X.g.train,2,sd),FUN="/")
Y.g.train = sweep(sweep(Y.g.train,2,colMeans(Y.g.train)),2,apply(Y.g.train,2,sd),FUN="/")

theta1.g = seq(0.2,0.4,.05)
theta2.g = 0.2
S.g = matrix(1,dim(Y.g.train)[1],dim(Y.g.train)[2])
W0.g = matrix(1,dim(X.g.train)[2],dim(Y.g.train)[2])
alpha.g = 0.01

rMSE.g = matrix(NA,length(theta1.g),4)
nMSE.g = rep(NA,length(theta1.g))
wR.g = rep(0,length(theta1.g))
for (a in 1:length(theta1.g))
{
  coef.g = gradientDescent(W0.g,alpha.g,theta1.g[a],theta2.g,S.g,X.g.train,Y.g.train,.1)
  X.g.test = X.interpolated[-train.indices,,1] 
  Y.g.test = Y.interpolated[-train.indices,1:4]
  Y.g.test = sweep(sweep(Y.g.test,2,colMeans(Y.g.test)),2,apply(Y.g.test,2,sd),FUN="/")
  X.g.test = sweep(sweep(X.g.test,2,colMeans(X.g.test)),2,apply(X.g.test,2,sd),FUN="/")
  predicted.g = X.g.test %*% coef.g
  
  for(i in 1:4)
  {
    rMSE.g[a,i] = sqrt(sum((Y.g.test[,i]-predicted.g[,i])^2)/dim(Y.g.test)[1])
  }
  nMSE.g[a] = sum((Y.g.test - predicted.g)^2)/(dim(Y.g.test)[1]*dim(Y.g.test)[2])
  for (i in 1:4)
  {
    wR.g[a] = wR.g[a] + cor(Y.g.test[,i],predicted.g[,i])*dim(Y.g.test)[1]
  }
  wR.g[a] = wR.g[a] / (dim(Y.g.test)[1] * dim(Y.g.test)[2])
}

