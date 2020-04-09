#testing - trying to see if this push works out

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#read in bootstrapped values for dose-response
exactbp<-read.csv('Exact_BetaPoisson_Bootstrap.csv')

require('gsl')

looplength<-10000

conc.1<-rep(NA,looplength)
conc.2<-rep(NA,looplength)
dose.1<-rep(NA,looplength)
dose.2<-rep(NA,looplength)
infect.1<-rep(NA,looplength)
infect.2<-rep(NA,looplength)
inhalation<-rep(NA,looplength)

  #1-4 micron (RNA/m^3)
  conc.1<-runif(looplength,min=916,max=1384)
  
  #>4 micron (RNA/m^3)
  conc.2<-runif(looplength,min=927,max=2000)
  
  #inhalation (original values in m^3/day... converted to per minute)
  inhalation<-runif(looplength,min=5.92,max=28.81)/(24*30*30) #minimum= (5th percentile), max (99th percentile)
  #for men and women, normal and obese, Table 6-6
  reduce<-sample(c(.25,.50,.75,.90,.95,.99),looplength,replace=TRUE)
  duration<-sample(c(30,5*60,15*60,1*60*60),looplength,replace=TRUE)
  dose.1<-conc.1*inhalation*duration*(1-reduce)
  dose.2<-conc.2*inhalation*duration*(1-reduce)
  
  for (i in 1:looplength){
    pair<-sample(c(1:length(exactbp$ln.alpha.)),1)
    infect.1[i]<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose.1[i], give=FALSE, strict=TRUE)
    infect.2[i]<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose.2[i], give=FALSE, strict=TRUE)
  }
  
  require(ggplot2)
  require(ggpubr)
  
  all<-data.frame(
    conc=c(conc.1,conc.2),
    reduce=rep(reduce,2),
    dose=c(dose.1,dose.2),
    infect=c(infect.1,infect.2),
    inhalation=rep(inhalation,2),
    type=c(rep("1-4 micron",looplength),rep(">4 micron",looplength)),
    duration=as.character(rep(duration,2))
  )
  
  all$duration[all$duration=="30"]<-"30 seconds"
  all$duration[all$duration=="300"]<-"5 minutes"
  all$duration[all$duration=="900"]<-"15 minutes"
  all$duration[all$duration=="3600"]<-"1 hour"
  
 A<-ggplot(data=all[all$type==">4 micron",],aes(x=conc,y=infect))+geom_point(aes(colour=as.character(reduce)),alpha=0.3)+
    facet_wrap(~duration)+
    scale_x_continuous(name=expression("Concentration (viral particles/m"^3*")"))+
    scale_y_continuous(trans="log10",name="Infection Risk")+
    scale_colour_discrete(name="Fraction Reduction")+
    theme_pubr()+theme(axis.text.x = element_text(angle = 90))+
    ggtitle("> 4 Micron")
 
B<-ggplot(data=all[all$type=="1-4 micron",],aes(x=conc,y=infect))+geom_point(aes(colour=as.character(reduce)),alpha=0.3)+
   facet_wrap(~duration)+
   scale_x_continuous(name=expression("Concentration (viral particles/m"^3*")"))+
   scale_y_continuous(trans="log10",name="Infection Risk")+
   scale_colour_discrete(name="Fraction Reduction")+
   theme_pubr()+theme(axis.text.x = element_text(angle = 90))+
   ggtitle("1-4 Micron")

windows()
ggarrange(A,B,common.legend=TRUE)
