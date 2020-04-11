
#functio requires definition of material type and exposure duration in minutes
COVIDmask<-function(material=c("100% cotton","scarf","tea towel","pillowcase","antimicrobial pillowcase",
                               "surgical mask","vacuum cleaner bag","cotton mix","linen","silk"),exposureduration){
  
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  #read in bootstrapped values for dose-response
  exactbp<-read.csv('Exact_BetaPoisson_Bootstrap.csv')
  
  require(gsl)
  require(truncdist)
  
  looplength<-10000
  
  conc.1<-rep(NA,looplength)
  conc.2<-rep(NA,looplength)
  dose.1<-rep(NA,looplength)
  dose.2<-rep(NA,looplength)
  infect<-rep(NA,looplength)
  inhalation<-rep(NA,looplength)
  
  #1-4 micron (RNA/m^3)
  conc.1<-runif(looplength,min=916,max=1384)
  
  #>4 micron (RNA/m^3)
  conc.2<-runif(looplength,min=927,max=2000)
  
  #inhalation (original values in m^3/day... converted to per minute)
  inhalation<-runif(looplength,min=5.92,max=28.81)/(24*30) #minimum= (5th percentile), max (99th percentile)
  #for men and women, normal and obese, Table 6-6
  
  #efficacies from Davies et al. (2013) MS2 values
  
  if (material=="100% cotton T-shirt"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.5085,sd=.1681)
  
  }else if (material=="scarf"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.4887,sd=.1977)
    
  }else if (material=="tea towel"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.7246,sd=.2260)
    
  }else if (material=="pillowcase"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.5713,sd=.1055)
    
  }else if (material=="antimicrobial pillowcase"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.6562,sd=.0764)
    
  }else if (material=="surgical mask"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.6890,sd=.0744)
    
  }else if (material=="vacuum cleaner bag"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.8952,sd=.0265)
    
  }else if (material=="cotton mix"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.7024,sd=.0008)
    
  }else if (material=="linen"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.6167,sd=.0241)
    
  }else{ #silk
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.5432,sd=.2949)
  }
  
  dose.1<-conc.1*inhalation*exposureduration*(1-reduce)
  dose.2<-conc.2*inhalation*exposureduration*(1-reduce)
  
  dose<-dose.1+dose.2
  for (i in 1:looplength){
    pair<-sample(c(1:length(exactbp$ln.alpha.)),1)
    infect[i]<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose[i], give=FALSE, strict=TRUE)
  }
  
  all<<-data.frame(
    conc=c(conc.1+conc.2),
    reduce=reduce,
    dose=dose,
    infect=infect,
    inhalation=inhalation,
    duration=as.character(rep(exposureduration,looplength)),
    materialtype=material
  )
  
}

#--------------- running sim and plotting -----------------------------------------
require(ggplot2)
require(ggpubr)

#scarf - 30 seconds and 15 minute exposure scenarios
COVIDmask(material="scarf",exposureduration=.5)
all.scarf.05<-all
COVIDmask(material="scarf",exposureduration=15)
all.scarf.15<-all

#linen - 30 seconds and 15 minute exposure scenarios
COVIDmask(material="linen",exposureduration=.5)
all.linen.05<-all
COVIDmask(material="linen",exposureduration=15)
all.linen.15<-all

#bind all scenarios into single data frame
all.materials<-rbind(all.scarf.05,all.scarf.15,all.linen.05,all.linen.15)

#plot
ggplot(all.materials)+
  geom_violin(aes(x=materialtype,y=infect,fill=materialtype),colour="black",draw_quantiles = c(.25,.5,.75))+
  #geom_jitter(aes(x=materialtype,y=infect),width=.1,alpha=.3)+
  theme_pubr()+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Infection Risk")+
  theme(legend.position = "none")+
  facet_wrap(~duration)
