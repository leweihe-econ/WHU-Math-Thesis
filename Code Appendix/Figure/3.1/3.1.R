
#Figure3.1

load(paste0(7,".rdata"))

library(scales)

p1<-hist(pi,breaks=seq(from=0,to=1,by=0.05))
p2<-hist(pscore,breaks=seq(from=0,to=1,by=0.05))
plot( p1, col='skyblue',xlim=c(0,1),border=F,
      main = 'Distribution of Propensity Scores',
      xlab = 'Propensity Scores')  # first histogram
plot( p2,col=scales::alpha('red',.45), xlim=c(0,1), add=T,border=F) 

legend("topright",                                    
       c("Correct Pi","Incorrect Pscore"),        
       fill=c("skyblue",scales::alpha('red',.45)),
       bty = 'n',border = NA) 

t.test(pi, pscore,
       alternative = c("two.sided"),
      paired = TRUE)