# Model PRS etc using relative risk
# Copyright 2019 David Curtis

# Update from code used in this paper:
# https://www.preprints.org/manuscript/201810.0171/v1

# This version models the PRS by assuming that the PRS has a linear effect on log of relative risk, close to logistic regression model
# For this, set expBeta which is the increase in relative risk associated with 1 unit increase in PRS
# Throughout, I have pretended that the odds ratio is the same as the risk ratio and refer to OR rather than RR
# I reckon this is good enough for a relatively rare trait and much easier to work with

library(ggplot2)

# Set these variables to the values you want

nSub=100000 # number of subjects to simulate
K=0.01 # prevalence
expBeta=1.7 # relative risk increase associated with 1 SD increase in PRS

# multiplot from Cookbook for R - http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

set.seed(1)

beta=log(expBeta)
risks=data.frame(matrix(0, ncol = 9, nrow =801))
colnames(risks)=c("PRS","Density","caseProportion","Risk","fractCont","fractCase","Specificity","Sensitivity","PPV")
risks$PRS=(-400:400)/100
risks$Density=dnorm(risks$PRS)
risks$caseProportion=exp(risks$PRS*beta)*risks$Density*0.01 # 0.01 is width of each PRS band
expAlpha=K/sum(risks$caseProportion)
risks$Risk=exp(risks$PRS*beta)*expAlpha
risks$fractCont=risks$Density*0.01*(1-risks$Risk)/(1-K)
risks$fractCase=risks$Density*0.01*risks$Risk/K
t=sprintf("Density distribution of PRS in cases and controls,\nK=%.4f, expBeta = %.3f",K,expBeta)
pDistRisk=ggplot(risks, aes(PRS)) +
  geom_line(aes(y=fractCont,colour="Controls")) + 
  geom_line(aes(y=fractCase,colour="Cases")) + 
  ggtitle(t) + 
  scale_x_continuous(breaks = seq(-4,4,by =1))  + 
  coord_cartesian(xlim=c(-4,4)) +
  scale_y_continuous(breaks = seq(0,1,by =0.1)) +
  theme_bw()

centiles=data.frame(matrix(0, ncol = 10, nrow =101))
colnames(centiles)=c("centile","PRS","controls","cases","Prevalence","Specificity","Sensitivity","PPV","AUC","OR")
centiles$centile=(0:100)
hundredth=nSub/100
totControl=nSub*(1-K)
totCase=nSub*K
totPercent=0
lastr=1
nControl=0
nCase=0
thisPercent=0
rr=2
for (r in 2:801)
{
  nControl=nControl+risks$fractCont[r]*totControl
  nCase=nCase+risks$fractCase[r]*totCase
  thisPercent=thisPercent+risks$Density[r]*0.01
  if (nControl+nCase>=hundredth || r==801) {
   thisLot=nControl+nCase
   centiles$controls[rr]=nControl*hundredth/thisLot
   centiles$cases[rr]=nCase*hundredth/thisLot
   centiles$PRS[rr]=(risks$PRS[lastr]+risks$PRS[r])/2
   lastr=r
   rr=rr+1
   totPercent=totPercent+0.01
   thisPercent=thisPercent-0.01
   nControl=nControl*(1-hundredth/thisLot)
   nCase=nCase*(1-hundredth/thisLot)
  }
}

centiles$cases[1]=0
centiles$controls[1]=hundredth
centiles$PRS[1]=centiles$PRS[2]
controlsSoFar=0
casesSoFar=0
centiles$Prevalence=centiles$cases/hundredth

totAUC=0
tCase=0.25
tControl=0.25
tOR=3
for (r in 2:101)
{
  controlsSoFar=controlsSoFar+centiles$controls[r]
  casesSoFar=casesSoFar+centiles$cases[r]
  centiles$Specificity[r]=controlsSoFar/totControl
  centiles$Sensitivity[r]=1-casesSoFar/totCase
  if ((1-centiles$Sensitivity[r])>tCase) {
    m=sprintf("Quartile for cases finishes at %d th centile",r)
    print(m)
    tCase=tCase+0.25
  }
  if (centiles$Specificity[r]>tControl) {
    m=sprintf("Quartile for controls finishes at %d th centile",r)
    print(m)
    tControl=tControl+0.25
  }

    AUC=(centiles$Sensitivity[r-1]+centiles$Sensitivity[r])*0.5*(centiles$Specificity[r]-centiles$Specificity[r-1])
    centiles$AUC[r]=AUC
    totAUC=totAUC+AUC
    centiles$PPV[r]=(totCase-casesSoFar)/((totCase-casesSoFar)+(totControl-controlsSoFar))
    if (r<101) {
      centiles$OR[r]=(totCase-casesSoFar)/(totControl-controlsSoFar)/(casesSoFar/controlsSoFar)
      if (centiles$OR[r]>tOR && r>50) {
        m=sprintf("OR reaches %d at %d th centile",tOR,r)
        print(m)
        tOR=tOR+1
      }
  }
}
centiles$AUC[1]=0
centiles$Specificity[1]=0
centiles$Sensitivity[1]=1
centiles$PPV[101]=centiles$PPV[100]

deciles=data.frame(matrix(0, ncol = 5, nrow =10))
colnames(deciles)=c("decile","PRS","controls","cases","relOR")
deciles$decile=(1:10)
for (r in 1:10)
{
  deciles$PRS[r]=centiles$PRS[(r-1)*10+6]
  deciles$controls[r]=0
  deciles$cases[r]=0
  for (rr in 1:10)
  {
    deciles$controls[r]=deciles$controls[r]+centiles$controls[(r-1)*10+1+rr]
    deciles$cases[r]=deciles$cases[r]+centiles$cases[(r-1)*10+1+rr]
  }
  if (r==1) {
    baseOdds=deciles$cases[r]/deciles$controls[r]
  }
  deciles$relOR[r]=(deciles$cases[r]/deciles$controls[r])/baseOdds
}

centiles=centiles[-1,]
t=sprintf("PRS decile relative ORs")
decORRisk=ggplot(deciles,aes(decile,relOR)) + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(1,10,by=1))  + expand_limits( y = 0) + scale_y_continuous(breaks=seq(0,100,1)) +
  ggtitle(t)
t=sprintf("OR dischotomised by PRS")
ORRisk=ggplot(centiles[-100,],aes(PRS,OR)) + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(-3,3,by=1))  + expand_limits( y = 0) + scale_y_continuous(breaks=seq(0,100,1)) +
  coord_cartesian(xlim=c(-4,4)) +
  ggtitle(t)
t=sprintf("Prevalence by centile of PRS")
PrevalenceRisk=ggplot(centiles,aes(centile,Prevalence)) + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(0,100,by =10)) + 
  scale_y_continuous(breaks = seq(0,K*5,by =K*0.2)) +
  ggtitle(t)
t=sprintf("Sensitivity, specificity and positive predictive value\nfor different values of PRS")
SpecSensRisk=ggplot(centiles,aes(PRS)) + 
  geom_line(aes(y=Specificity,colour="Specificity")) + 
  geom_line(aes(y=Sensitivity,colour="Sensitivity")) + 
  geom_line(aes(y=PPV,colour="PPV")) + 
  ggtitle(t) + 
  scale_x_continuous(breaks = seq(-3,3,by =1))  + 
  coord_cartesian(xlim=c(-4,4)) +
  scale_y_continuous(breaks = seq(0,1,by =0.1)) +
  ylab("Sensitivity, specificity, PPV") +
  theme_bw()
t=sprintf("Sensitivity, specificity and positive predictive value\nfor different centiles of PRS")
SpecSensCentRisk=ggplot(centiles,aes(centile)) + 
  geom_line(aes(y=Specificity,colour="Specificity")) + 
  geom_line(aes(y=Sensitivity,colour="Sensitivity")) + 
  geom_line(aes(y=PPV,colour="PPV")) + 
  ggtitle(t) + scale_x_continuous(breaks = seq(-3,3,by =1))  + 
  scale_x_continuous(breaks = seq(0,100,by =10)) + 
  scale_y_continuous(breaks = seq(0,1,by =0.1)) +
  ylab("Sensitivity, specificity, PPV") +
  theme_bw()
t=sprintf("Sensitivity against specificity,\nAUC = %f",totAUC)
ROCRisk=ggplot(centiles,aes(Specificity,Sensitivity)) + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(0,1,by =0.1)) + 
  scale_y_continuous(breaks = seq(0,1,by =0.1)) +
  ggtitle(t)
r <- vector("list", 6)
r[[1]]=pDistRisk
r[[2]]=decORRisk
r[[3]]=PrevalenceRisk
r[[4]]=ORRisk
r[[5]]=SpecSensCentRisk
r[[6]]=ROCRisk
dev.new()
multiplot(plotlist=r,cols=3)

