#Code for estimating infection risks for different mask materials

#clear environemnt
rm(list = ls())

#set seed
set.seed(34)

#set working directory to location of this and other files
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read in bootstrapped values for dose-response
exactbp<-read.csv('Exact_BetaPoisson_Bootstrap.csv')

#function requires definition of material type and exposure duration in minutes
COVIDmask<-function(material=c("none","100% cotton","scarf","tea towel","pillowcase","antimicrobial pillowcase",
                               "surgical mask","vacuum cleaner bag","cotton mix","linen","silk"),
                    exposureduration,
                    RNAinfective){
  
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
  conc.1<-runif(looplength,min=916,max=1384)*RNAinfective
  
  #>4 micron (RNA/m^3)
  conc.2<-runif(looplength,min=927,max=2000)*RNAinfective
  
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
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.6890,sd=.0744)
    
  }else if (material=="surgical mask"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.8952,sd=.0265)
    
  }else if (material=="vacuum cleaner bag"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.8595,sd=.0155)
    
  }else if (material=="cotton mix"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.7024,sd=.0008)
    
  }else if (material=="linen"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.6167,sd=.0241)
    
  }else if (material=="silk"){
    reduce<-rtrunc(looplength,"norm",a=0,b=1,mean=.5432,sd=.2949)
  }else{
    reduce<-rep(0,looplength)
  }
  
  dose.1<-conc.1*inhalation*exposureduration*(1-reduce)
  dose.2<-conc.2*inhalation*exposureduration*(1-reduce)
  
  dose<-dose.1+dose.2
  for (i in 1:looplength){
    pair<-sample(c(1:length(exactbp$ln.alpha.)),1)
    infect[i]<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose[i], give=FALSE, strict=TRUE)
  }
  
  all.param<<-data.frame(
    conc=c(conc.1+conc.2),
    reduce=reduce,
    dose=dose,
    infect=infect,
    inhalation=inhalation,
    duration=as.character(rep(exposureduration,looplength)),
    materialtype=material,
    RNAinfect=rep(RNAinfective,looplength)
  )

}

#--------------- function for running simulation scenarios per material type -----------------------------------------

durationandRNA<-function(exposureduration,RNAinfective){
  #scarf - 30 seconds and 15 minute exposure scenarios
  COVIDmask(material="scarf",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.scarf<-all.param

  #linen - 30 seconds and 15 minute exposure scenarios
  COVIDmask(material="linen",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.linen<-all.param
  
  #t shirt
  COVIDmask(material="100% cotton T-shirt",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.tshirt<-all.param
  
  #cotton mix
  COVIDmask(material="cotton mix",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.cottonmix<-all.param
  
  #antimicrobial pillowcase
  COVIDmask(material="antimicrobial pillowcase",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.antimicrobepillowcase<-all.param
  
  #pillowcase
  COVIDmask(material="pillowcase",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.pillowcase<-all.param

  #vacuum cleaner bag
  COVIDmask(material="vacuum cleaner bag",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.vacuum<-all.param
  
  #surgical mask
  COVIDmask(material="surgical mask",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.surgicalmask<-all.param

  #tea towel
  COVIDmask(material="tea towel",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.teatowel<-all.param
  
  #silk
  COVIDmask(material="silk",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.silk<-all.param
  
  #none
  COVIDmask(material="none",exposureduration=exposureduration,RNAinfective=RNAinfective)
  all.none<-all.param
  
  #bind all scenarios into single data frame
  all.materials<<-rbind(all.scarf,
                       all.linen,
                       all.tshirt,
                       all.cottonmix,
                       all.antimicrobepillowcase,
                       all.pillowcase,
                       all.vacuum,
                       all.surgicalmask,
                       all.teatowel,
                       all.silk,
                       all.none
  )
  
}


#running simulations for different material types, exposure times, and % of RNA assumed to be infective-----------------------------------------
exposuretimes<-c(.5,20)
RNAinfect<-c(.001,.01,.1)
for (i in 1:length(exposuretimes)){
  for (j in 1:length(RNAinfect)){
    durationandRNA(exposureduration=exposuretimes[i],RNAinfective=RNAinfect[j])
    if (j==1 & i==1){
      all.materials.total<-all.materials 
    }else{
      all.materials.total<-rbind(all.materials,all.materials.total)
    }
  }
}


#Update for labeling purposes
all.materials.total$duration<-as.character(all.materials.total$duration)
all.materials.total$duration[all.materials.total$duration==.5]<-"30 Seconds"
all.materials.total$duration[all.materials.total$duration==20]<-"20 Minutes"

all.materials.total$RNAinfect[all.materials.total$RNAinfect==.001]<-"0.1% Infective"
all.materials.total$RNAinfect[all.materials.total$RNAinfect==.01]<-"1% Infective"
all.materials.total$RNAinfect[all.materials.total$RNAinfect==.1]<-"10% Infective"


#figure 1 --------------------------------------------------------------------------------------------------------------------------------
require(ggplot2)
require(ggpubr)

windows()
ggplot(all.materials.total)+
  geom_violin(aes(x=materialtype,y=infect,fill=materialtype),colour="black",draw_quantiles = c(.25,.5,.75))+
  theme_pubr()+
  scale_x_discrete(name="")+
  scale_y_continuous(name="Infection Risk",trans="log10")+
  theme(legend.position = "none")+
  coord_flip()+
  facet_wrap(RNAinfect~duration,ncol=2)

#summary statistics ----------------------------------------------------------------

materials<-c("none","silk","tea towel","surgical mask","vacuum cleaner bag","pillowcase",
             "antimicrobial pillowcase", "cotton mix", "100% cotton T-shirt","linen","scarf")

none.mean.20min<-mean(all.materials.total$infect[all.materials.total$materialtype=="none" & all.materials.total$duration=="20 Minutes"])
none.mean.30sec<-mean(all.materials.total$infect[all.materials.total$materialtype=="none" & all.materials.total$duration=="30 Seconds"])

matrix.reduce<-matrix(ncol=length(materials),nrow=length(exposuretimes))
colnames(matrix.reduce)<-materials
rownames(matrix.reduce)<-c("20 Minutes","30 Seconds")

for (i in 1:length(materials)){
  material.20min<-mean(all.materials.total$infect[all.materials.total$materialtype==materials[i] & all.materials.total$duration=="20 Minutes"])
  material.30sec<-mean(all.materials.total$infect[all.materials.total$materialtype==materials[i] & all.materials.total$duration=="30 Seconds"])
  
  matrix.reduce[1,i]<-(none.mean.20min-material.20min)/none.mean.20min*100
  matrix.reduce[2,i]<-(none.mean.30sec-material.30sec)/none.mean.30sec*100
}

View(matrix.reduce)
