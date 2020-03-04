
library(Matching)

# Mahalanobis Distance

ATTmh<-c()

for (r in 1:1000) {
  
  load(paste0(r,".rdata"))
  attach(dataset)
  
  X<-dataset[,c('age','ed','black','hisp','married'
                ,'nodeg','re74','re75','u74','u75')]
  
  res<-Match(Y=dataset$y,Tr=dataset$treated, X=X,estimand = "ATT",
             Weight = 2, weights=rep(1,445), ties = FALSE)
  
  ATTmh[r]<-res$est.noadj
}



#Propensity Scores

ATTep<-c()

for (r in 1:1000) {
  
  load(paste0(r,".rdata"))
  attach(dataset)
  
  X<-dataset[,c('age','ed','black','hisp','married'
                ,'nodeg','re74','re75','u74','u75')]
  
  res<-Match(Y=dataset$y,Tr=dataset$treated, X=dataset$pscore,
             estimand = "ATT",
             Weight = 2)
  
  ATTep[r]<-res$est
}


