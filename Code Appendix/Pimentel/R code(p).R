setwd("~/Wuhan University/WHU thesis/math/pscore")

ATTp<-c()

for ( r in 1:1000) {

load(paste0(r,".rdata"))

library(Matching)
library(rgenoud)

library(rcbalance)
library('DOS')
library(optmatch)

data1<-dataset[dataset$pscore>1/3,]
data2<-dataset[dataset$pscore<=1/3&dataset$pscore>1/4,]
data3<-dataset[dataset$pscore<=1/4,]


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
  
  distmat1<-build.dist.struct(data1$treated, data1[,2:9], 
                            calip.option = "propensity",caliper = 100)
  
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
  
  
  try(distmat2<-build.dist.struct(data2$treated, data2[,2:9], 
                              calip.option = "propensity",caliper = 100), silent = TRUE)
  
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
  
  distmat3<-build.dist.struct(data3$treated, data3[,2:9], 
                              calip.option = "propensity",caliper = 100)
  
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
  
  treat1<-data1[data1$treated==1,c('age','ed','black','hisp','married'
                                   ,'nodeg','re74','re75','u74','u75')]
  control1<-data1[data1$treated==0,c('age','ed','black','hisp','married'
                                     ,'nodeg','re74','re75','u74','u75')]
  
  no.control1<-max(c(1,min(c( 2, 5,
                             floor(nrow(dmat1)/ncol(dmat1))))))
  
  
  if(nrow(dmat1)/ncol(dmat1)<1){
    match1<-rcbalance(distmat1,near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat1, control.info = control1,
                      tol = 0.1, k=no.control1, exclude.treated=TRUE)
  }else{
    match1<-rcbalance(distmat1,near.exact = c('age','ed','black','hisp','married'
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
  
  treat2<-data2[data2$treated==1,c('age','ed','black','hisp','married'
                                   ,'nodeg','re74','re75','u74','u75')]
  control2<-data2[data2$treated==0,c('age','ed','black','hisp','married'
                                     ,'nodeg','re74','re75','u74','u75')]
  
  no.control2<-max(c(1,min(c( 3, 5,
                              floor(nrow(dmat2)/ncol(dmat2))))))
  
  
  if(nrow(dmat2)/ncol(dmat2)<1){
    try(match2<-rcbalance(distmat2,near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat2, control.info = control2,
                      tol = 0.1, k=no.control2, exclude.treated=TRUE),silent = TRUE)
  }else{
    try(match2<-rcbalance(distmat2,near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat2, control.info = control2,
                      tol = 0.1, k=no.control2), silent = TRUE)
  }
  
  try(matches2<-match2$matches,silent = TRUE)
  
  
  #who is matched with who
  
  try(matchmat2<-matrix(0, nrow(matches2), 
                    ncol =no.control2+1),silent = TRUE)
  try(matchmat2[,1]<-as.integer(colnames(dmat2)[as.integer(rownames(matches2))]),silent = TRUE)
  
  
  try(for (i in 1:nrow(matchmat2)) {
    for (j in 2:(no.control2+1)) {
      matchmat2[i,j]<-as.integer(rownames(dmat2))[match2$matches[i,j-1]]
    }
  }, silent = TRUE)
  
}

#############
#Group3

if(length(data3$treated[data3$treated==0])!=0&
   length(data3$treated[data3$treated==1])!=0){
  
  
  #matchvec=pairmatch(t(dmat1),controls=2)
  #matchvec
  
  treat3<-data3[data3$treated==1,c('age','ed','black','hisp','married'
                                   ,'nodeg','re74','re75','u74','u75')]
  control3<-data3[data3$treated==0,c('age','ed','black','hisp','married'
                                     ,'nodeg','re74','re75','u74','u75')]
  
  no.control3<-max(c(1,min(c( 4, 5,
                              floor(nrow(dmat3)/ncol(dmat3))))))
  
  
  if(nrow(dmat3)/ncol(dmat3)<1){
    match3<-rcbalance(distmat3,near.exact = c('age','ed','black','hisp','married'
                                              ,'nodeg','re74','re75','u74','u75'),
                      fb.list = c('married'),treated.info = treat3, control.info = control3,
                      tol = 0.1, k=no.control3, exclude.treated=TRUE)
  }else{
    match3<-rcbalance(distmat3,near.exact = c('age','ed','black','hisp','married'
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
   length(data1$treated[data1$treated==1])!=0&
   exists('matchmat3')){
  
  ytreat1<-dataset$y[matchmat1[,1]]
  ycontrol1<-c()
  
  for (i in 1:length(ytreat1)) {
    ycontrol1[i]<-mean(dataset$y[matchmat1[i,2:(no.control1+1)]])
  }
  
}else{
  ytreat1<-c()
  ycontrol1<-c()
}

#Group2

if(length(data2$treated[data2$treated==0])!=0&
   length(data2$treated[data2$treated==1])!=0&
   exists('matchmat2')){
  
  ytreat2<-dataset$y[matchmat2[,1]]
  ycontrol2<-c()
  
  for (i in 1:length(ytreat2)) {
    ycontrol2[i]<-mean(dataset$y[matchmat2[i,2:(no.control2+1)]])
  }
  
}else{
  ytreat2<-c()
  ycontrol2<-c()
}

#Group3

if(length(data3$treated[data3$treated==0])!=0&
   length(data3$treated[data3$treated==1])!=0&
   exists('matchmat3')){
  
  ytreat3<-dataset$y[matchmat3[,1]]
  ycontrol3<-c()
  
  for (i in 1:length(ytreat3)) {
    ycontrol3[i]<-mean(dataset$y[matchmat3[i,2:(no.control3+1)]])
  }
  
}else{
  ytreat3<-c()
  ycontrol3<-c()
}

ytreat<-c(ytreat1,ytreat2,ytreat3)
ycontrol<-c(ycontrol1,ycontrol2,ycontrol3)

TT<-ytreat-ycontrol

ATTp[r]<-mean(TT)

}

# replace NA
ATTp[is.na(ATTp)]=mean(ATTp,na.rm=T)

