
#################################################
# Data Read #
#################################################

library(Matrix)
library(haven)
setwd("~/Wuhan University/WHU thesis/math/data")
dataset <- read_stata('nswre74.dta')
View(dataset)

#################################################
# Generate Outcome #
#################################################

attach(dataset)

dataset$index<-seq(1,445)

dataset$cy<-1000*treat+0.1*exp(0.7*log(re74+0.01)+0.7*log(re75+0.01))

dataset$u74<-ifelse(re74==0,1,0)
dataset$u75<-ifelse(re75==0,1,0)
  
dataset$mu<-1+1.428*10^(-4)*age2-2.918*10^(-3)*ed^2-0.2275*black-0.8276*hisp
            +0.2071*married-0.8232*nodeg-1.236*10^(-9)*re74^2+5.865*10^(-10)*re75^2
            -0.4328*dataset$u74-0.3804*dataset$u75


dataset$pi<-1/(1+
                exp(-(20+0.5*dataset$mu+0.01*age2-0.3*ed^2
                -0.01*log(re74+0.01)^2+0.01*log(re75+0.01)^2-2*nodeg+1.5*married)))

summary(dataset$pi)


#################################################
# Simulate 1000 datasets #
#################################################
setwd("~/Wuhan University/WHU thesis/math/data copy")

for ( r in 1:1000) {
  
  set.seed(r)
  
  for (i in 1:445) {
    dataset$treated[i]<-rbinom(1,1,dataset$pi[i])
  }
  
  propscore.model=glm(dataset$treated~age+ed+black+hisp+married
                      +nodeg+re74+re75+u74+u75,family=binomial,
                      x=TRUE,y=TRUE,data=dataset)
  
  dataset$pscore<-1/(1+exp(-predict(propscore.model)))
  
  set.seed(r+1)
  
  dataset$epsilon<-10*rnorm(length(re75))
  
  dataset$y<-dataset$cy+dataset$epsilon
  
  save(dataset,file = paste0(r,".rdata"))
  
}

#################################################
# Weighted Matrix #
#################################################

#install.packages("rgenoud")
#install.packages("Matching", dependencies=TRUE)


library(Matching)
library(rgenoud)

library(rcbalance)
library('DOS')
library(optmatch)

library(flexclust)

flexclust::distCor
function (x, centers) 
{
  z <- matrix(0, nrow(x), ncol = nrow(centers))
  for (i in 1:nrow(x)) {
    for (j in 1:nrow(centers)) {
      z[i,j]<-(t(x[i,]-centers[j,])%*%t(solve(Cholesky(S)))%*%
                 W%*%solve(Cholesky(S))%*%(x[i,]-centers[j,]))^0.5
    }
  }
  z
}

#ATT<-c()
#aATT<-c()

for (r in 1:1000) {
  
  load(paste0(r,".rdata"))
  attach(dataset)
  
  result<-GenMatch(treat, dataset[,c('age','ed','black','hisp','married'
                                     ,'nodeg','re74','re75','u74','u75')],
                   min.weight=0, max.weight=1000, max.generations=11, 
                   wait.generations=3,estimand="ATE")
  
  Weight<-result[3]
  W<-matrix(unlist(Weight),nrow=10)
  
  
  #mout <- Match(Y=y, Tr=treat, dataset[,c('age','ed','black','hisp','married'
  #                                        ,'nodeg','re74','re75','u74','u75')],
  #             estimand="ATT", Weight.matrix=result)
  #summary(mout)
  
  #################################################
  # Distance Measure #
  #################################################
  
  X<-dataset[,c('age','ed','black','hisp','married'
                ,'nodeg','re74','re75','u74','u75')]
  
  S<-cov(X)
  
  
  #################################################
  # Clustering #
  #################################################
  
  res = kcca(X, 3, family=kccaFamily(dist=distCor))
  attributes(res)
  
  dataset$group<-attributes(res)$cluster
  
  data1<-dataset[dataset$group==1,]
  data2<-dataset[dataset$group==2,]
  data3<-dataset[dataset$group==3,]
  
  
  #################################################
  # Distance Matrix #
  #################################################
  
  if(length(data1$treated[data1$treated==0])==0|
     length(data1$treated[data1$treated==1])==0)
    {
    ytreat1<-c()
    ycontrol1<-c()
  }else{
    # Group1
    
    dmat1<-matrix(0, length(data1$treated[data1$treated==0]), 
                  ncol =length(data1$treated[data1$treated==1]))
    
    colnames(dmat1)<-c(as.character(data1$index[data1$treated==1]))
    rownames(dmat1)<-c(as.character(data1$index[data1$treated==0]))
    
    
    for (i in 1:nrow(dmat1)) {
      for (j in 1:ncol(dmat1)) {
        dmat1[i,j]<-(as.matrix(X[as.numeric(rownames(dmat1)[i]),]-X[as.numeric(colnames(dmat1)[j]),])%*%
                       t(solve(chol(S)))%*%
                       W%*%solve(chol(S))%*%
                       t(as.matrix(X[as.numeric(rownames(dmat1)[i]),]-X[as.numeric(colnames(dmat1)[j]),])))^0.5
      }
    }
    
  }
  
  #Group2
  if(length(data2$treated[data2$treated==0])==0|
     length(data2$treated[data2$treated==1])==0)
  {
    ytreat2<-c()
    ycontrol2<-c()
  }else{
  
  dmat2<-matrix(0, length(data2$treated[data2$treated==0]), 
                ncol =length(data2$treated[data2$treated==1]))
  
  
  colnames(dmat2)<-c(as.character(data2$index[data2$treated==1]))
  rownames(dmat2)<-c(as.character(data2$index[data2$treated==0]))
  
  
  for (i in 1:nrow(dmat2)) {
    for (j in 1:ncol(dmat2)) {
      dmat2[i,j]<-(as.matrix(X[as.numeric(rownames(dmat2)[i]),]-X[as.numeric(colnames(dmat2)[j]),])%*%
                     t(solve(chol(S)))%*%
                     W%*%solve(chol(S))%*%
                     t(as.matrix(X[as.numeric(rownames(dmat2)[i]),]-X[as.numeric(colnames(dmat2)[j]),])))^0.5
      }
    }
  
  }
  
  
  
  #Group3
  if(length(data3$treated[data3$treated==0])==0|
     length(data3$treated[data3$treated==1])==0)
  {
    ytreat3<-c()
    ycontrol3<-c()
  }else{
  dmat3<-matrix(0, length(data3$treated[data3$treated==0]), 
                ncol =length(data3$treated[data3$treated==1]))
  
  
  colnames(dmat3)<-c(as.character(data3$index[data3$treated==1]))
  rownames(dmat3)<-c(as.character(data3$index[data3$treated==0]))
  
  
  for (i in 1:nrow(dmat3)) {
    for (j in 1:ncol(dmat3)) {
      dmat3[i,j]<-(as.matrix(X[as.numeric(rownames(dmat3)[i]),]-X[as.numeric(colnames(dmat3)[j]),])%*%
                     t(solve(chol(S)))%*%
                     W%*%solve(chol(S))%*%
                     t(as.matrix(X[as.numeric(rownames(dmat3)[i]),]-X[as.numeric(colnames(dmat3)[j]),])))^0.5
    }
  }
  
  }
  
  
  #################################################
  # Fine Balance Matching #
  #################################################
  
  
  
  #############
  #Group1
  
  if(length(data1$treated[data1$treated==0])!=0&
     length(data1$treated[data1$treated==1])!=0){
  
  #matchvec=pairmatch(t(dmat1),controls=2)
  #matchvec
  
  data1$hed<-ifelse(data1$ed>10,1,0)
  
  treat1<-data1[data1$treated==1,c('age','ed','black','hisp','married'
                                   ,'nodeg','re74','re75','u74','u75')]
  control1<-data1[data1$treated==0,c('age','ed','black','hisp','married'
                                     ,'nodeg','re74','re75','u74','u75')]
  
  no.control1<-max(c(1,min(c(floor(mean((1-data1$pscore)/data1$pscore)),2,
                             floor(nrow(dmat1)/ncol(dmat1))))))
  
  
  if(nrow(dmat1)/ncol(dmat1)<1){
    match1<-rcbalance(t(dmat1),near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat1, control.info = control1,
                      tol = 0.1, k=no.control1, exclude.treated=TRUE)
  }else{
    match1<-rcbalance(t(dmat1),near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat1, control.info = control1,
                      tol = 0.1, k=no.control1)
  }
  
  matches1<-match1$matches
  
  #who is matched with who
  
  matchmat1<-matrix(0, nrow(matches1), 
                    ncol =no.control1+1)
  matchmat1[,1]<-as.integer(colnames(dmat1)[as.integer(rownames(matches1))])
  
  
  for (i in 1:nrow(matchmat1)) {
    for (j in 2:(no.control1+1)) {
      matchmat1[i,j]<-as.integer(rownames(dmat1))[match1$matches[i,j-1]]
      }
    }
  }
  
  #############
  #Group2
  
  if(length(data2$treated[data2$treated==0])!=0&
     length(data2$treated[data2$treated==1])!=0){
    
  
  #matchvec=pairmatch(t(dmat1),controls=2)
  #matchvec
  
  data2$hed<-ifelse(data2$ed>10,1,0)
  
  treat2<-data2[data2$treated==1,c('age','ed','black','hisp','married'
                                   ,'nodeg','re74','re75','u74','u75')]
  control2<-data2[data2$treated==0,c('age','ed','black','hisp','married'
                                     ,'nodeg','re74','re75','u74','u75')]
  
  no.control2<-max(c(1,min(c(floor(mean((1-data2$pscore)/data2$pscore)),2,
                             floor(nrow(dmat2)/ncol(dmat2))))))
  
  if(nrow(dmat2)/ncol(dmat2)<1){
    match2<-rcbalance(t(dmat2),near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat2, control.info = control2,
                      tol = 0.1, k=no.control2, exclude.treated=TRUE)
  }else{
    match2<-rcbalance(t(dmat2),near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat2, control.info = control2,
                      tol = 0.1, k=no.control2)
  }
  
  matches2<-match2$matches
  
  #who is matched with who
  
  matchmat2<-matrix(0, nrow(matches2), 
                    ncol =no.control2+1)
  matchmat2[,1]<-as.integer(colnames(dmat2)[as.integer(rownames(matches2))])
  
  
  for (i in 1:nrow(matchmat2)) {
    for (j in 2:(no.control2+1)) {
      matchmat2[i,j]<-as.integer(rownames(dmat2))[match2$matches[i,j-1]]
    }
  }
  
  }
  
  #############
  #Group3
  
  if(length(data3$treated[data3$treated==0])!=0&
     length(data3$treated[data3$treated==1])!=0){
    
  
  #matchvec=pairmatch(t(dmat1),controls=2)
  #matchvec
  
  
  data3$hed<-ifelse(data3$ed>10,1,0)
  
  treat3<-data3[data3$treated==1,c('age','ed','black','hisp','married'
                                   ,'nodeg','re74','re75','u74','u75','hed')]
  control3<-data3[data3$treated==0,c('age','ed','black','hisp','married'
                                     ,'nodeg','re74','re75','u74','u75','hed')]
  
  no.control3<-max(c(1,min(c(floor(mean((1-data3$pscore)/data3$pscore)),2,
                             floor(nrow(dmat3)/ncol(dmat3))))))
  
  
  if(nrow(dmat3)/ncol(dmat3)<1){
    match3<-rcbalance(t(dmat3),near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat3, control.info = control3,
                      tol = 0.1, k=no.control3, exclude.treated=TRUE)
  }else{
    match3<-rcbalance(t(dmat3),near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat3, control.info = control3,
                      tol = 0.1, k=no.control3)
  }
  
  matches3<-match3$matches
  
  #who is matched with who
  
  matchmat3<-matrix(0, nrow(matches3), 
                    ncol =no.control3+1)
  matchmat3[,1]<-as.integer(colnames(dmat3)[as.integer(rownames(matches3))])
  
  
  for (i in 1:nrow(matchmat3)) {
    for (j in 2:(no.control3+1)) {
      matchmat3[i,j]<-as.integer(rownames(dmat3))[match3$matches[i,j-1]]
    }
  }
  
  }
  
  #################################################
  # Inference #
  #################################################
  
  #Group1
  
  if(length(data1$treated[data1$treated==0])!=0&
     length(data1$treated[data1$treated==1])!=0){
    
    ytreat1<-dataset$y[matchmat1[,1]]
    ycontrol1<-c()
    
    for (i in 1:length(ytreat1)) {
      ycontrol1[i]<-mean(dataset$y[matchmat1[i,2:(no.control1+1)]])
    }
    
  }
  
  #Group2
  
  if(length(data2$treated[data2$treated==0])!=0&
     length(data2$treated[data2$treated==1])!=0){
    
    ytreat2<-dataset$y[matchmat2[,1]]
    ycontrol2<-c()
    
    for (i in 1:length(ytreat2)) {
      ycontrol2[i]<-mean(dataset$y[matchmat2[i,2:(no.control2+1)]])
    }
    
  }
  
  #Group3
  
  if(length(data3$treated[data3$treated==0])!=0&
     length(data3$treated[data3$treated==1])!=0){
    
    ytreat3<-dataset$y[matchmat3[,1]]
    ycontrol3<-c()
    
    for (i in 1:length(ytreat3)) {
      ycontrol3[i]<-mean(dataset$y[matchmat3[i,2:(no.control3+1)]])
    }
    
  }
  
  ytreat<-c(ytreat1,ytreat2,ytreat3)
  ycontrol<-c(ycontrol1,ycontrol2,ycontrol3)

  TT<-ytreat-ycontrol
  
  ATT[r]<-mean(TT)
  
}


ATTnew<-ATT



