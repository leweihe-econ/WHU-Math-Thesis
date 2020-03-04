
#boxplot 3.2

bx<-c(ATTnew,ATTmh)
bx<-data.frame(bx)
bx$Method[1:1000]<-rep('New',1000)
bx$Method[1001:2000]<-rep('Mahalanobis',1000)

colnames(bx)<-c('ATT','Method')

library(ggplot2)

p <- ggplot(bx, aes(x=Method, y=ATT)) + 
  geom_boxplot() + 
  coord_flip()+
  labs(title="Comparision of Algorithms")+
  theme(plot.title = element_text(hjust = 0.5)) 

p

# Table 3.2

RMSEnew<-sqrt(sum((ATTnew-1000)^2)/1000)
RMSEp<-sqrt(sum((ATTp-1000)^2)/1000)
RMSEmh<-sqrt(sum((ATTmh-1000)^2)/1000)
RMSEep<-sqrt(sum((ATTep-1000)^2)/1000)

Biasnew<-mean(ATTnew)-1000
Biasp<-mean(ATTp)-1000
Biasmh<-mean(ATTmh)-1000
Biasep<-mean(ATTep)-1000

  