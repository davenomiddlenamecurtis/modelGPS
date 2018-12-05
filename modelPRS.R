# Model PRS etc contributing to variance in liability threshold model0
# Copyright 2018 David Curtis

library(ggplot2)

# Set these variables to the values you want

nSub=100000 # number of subjects to simulate
pVar=0.06 # proportion of variance explained
K=0.04 # prevalence

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
lT=qnorm(1-K)
sP=sqrt(pVar)
sR=sqrt(1-pVar)
subs=data.frame(matrix(0, ncol = 8, nrow =nSub))
colnames(subs)=c("rawP","PRS","Residual","Liability","Phenotype","Specificity","Sensitivity","PPV")
subs$rawP=rnorm(nSub,sd=sP)
subs$Residual=rnorm(nSub,sd=sR)
subs$Liability=subs$rawP+subs$Residual
subs$PRS=subs$rawP/sP
subs$Phenotype[subs$Liability<lT]="controls"
subs$Phenotype[subs$Liability>=lT]="cases"
subs=subs[order(subs$PRS),]
l <- vector("list", 6)
t=sprintf("Distribution of liability in cases and controls,\nK = %.4f, liability threshold = %.2f",K,lT)
lFreq=ggplot(subs, aes(x=Liability,color=Phenotype)) + geom_freqpoly(binwidth=0.05) + theme_bw() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  ggtitle(t) + scale_x_continuous(breaks = seq(-3,3,by =1)) + coord_cartesian(xlim=c(-4,4))
t=sprintf("Frequency distribution of PRS in cases and controls,\nK=%.4f, pVar = %.3f",K,pVar)
pFreq=ggplot(subs, aes(x=PRS,color=Phenotype)) + geom_freqpoly(binwidth=0.1) + theme_bw() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  ggtitle(t) + scale_x_continuous(breaks = seq(-3,3,by =1))  + coord_cartesian(xlim=c(-4,4))
t=sprintf("Density distribution of PRS in cases and controls,\nK=%.4f, pVar = %.3f",K,pVar)
t=sprintf("Density distribution of PRS in cases and controls,\nK=%.4f, pVar = %.3f",K,pVar)
pDist=ggplot(subs, aes(x=PRS,color=Phenotype)) + geom_density() + theme_bw() +
  ggtitle(t) + scale_x_continuous(breaks = seq(-3,3,by =1))  + coord_cartesian(xlim=c(-4,4))
# l[[1]]=lFreq
l[[1]]=pFreq
l[[2]]=pDist
centiles=data.frame(matrix(0, ncol = 10, nrow =101))
colnames(centiles)=c("centile","PRS","controls","cases","Prevalence","Specificity","Sensitivity","PPV","AUC","OR")
centiles$centile=(0:100)
hundredth=nSub/100
totControl=0
totCase=0
for (r in 2:101)
{
  centiles$PRS[r]=subs$PRS[nSub/100*(r-1)]
  part=subs$Phenotype[(hundredth*(centiles$centile[r]-1)+1):(hundredth*(centiles$centile[r]))]
  nControl=length(which(part=="controls"))
  centiles$controls[r]=nControl
  totControl=totControl+nControl
  nCase=length(which(part=="cases"))
  centiles$cases[r]=nCase
  totCase=totCase+nCase
}
centiles$cases[1]=0
centiles$controls[1]=0
centiles$PRS[1]=centiles$PRS[2]
controlsSoFar=0
casesSoFar=0
centiles$Prevalence=centiles$cases/hundredth
totAUC=0
tCase=0.25
tControl=0.25
tOR=3
for (r in 1:101)
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
  if (r==1) {
    centiles$AUC[r]=0 
  } else {
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
}
centiles$PPV[101]=centiles$PPV[100]
t=sprintf("OR against PRS")
OR=ggplot(centiles,aes(PRS,OR)) + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(-3,3,by =1))  + 
  coord_cartesian(xlim=c(-4,4)) +
  ggtitle(t)
t=sprintf("Prevalence by centile of PRS")
Prevalence=ggplot(centiles,aes(centile,Prevalence)) + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(0,100,by =10)) + 
  scale_y_continuous(breaks = seq(0,1,by =0.05)) +
  ggtitle(t)
t=sprintf("Sensitivity, specificity and positive predictive value\nfor different values of PRS")
SpecSens=ggplot(centiles,aes(PRS)) + 
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
SpecSensCent=ggplot(centiles,aes(centile)) + 
  geom_line(aes(y=Specificity,colour="Specificity")) + 
  geom_line(aes(y=Sensitivity,colour="Sensitivity")) + 
  geom_line(aes(y=PPV,colour="PPV")) + 
  ggtitle(t) + scale_x_continuous(breaks = seq(-3,3,by =1))  + 
  scale_x_continuous(breaks = seq(0,100,by =10)) + 
  scale_y_continuous(breaks = seq(0,1,by =0.1)) +
  ylab("Sensitivity, specificity, PPV") +
  theme_bw()
t=sprintf("Sensitivity against specificity,\nAUC = %f",totAUC)
ROC=ggplot(centiles,aes(Specificity,Sensitivity)) + geom_line() + theme_bw() +
  scale_x_continuous(breaks = seq(0,1,by =0.1)) + 
  scale_y_continuous(breaks = seq(0,1,by =0.1)) +
  ggtitle(t)
l[[3]]=Prevalence
l[[4]]=OR
l[[5]]=SpecSensCent
l[[6]]=ROC
multiplot(plotlist=l,cols=3)


