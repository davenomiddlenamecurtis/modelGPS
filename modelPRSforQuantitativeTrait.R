# Plot quantitative trait against PRS given correlation and inter-decile spread
# Set these variables to the values you want
# These values based on Figure 2B in https://www.cell.com/cell/pdf/S0092-8674(19)30290-9.pdf

r=0.29 # correlation coefficient
pVar=r*r # proportion of variance explained
lowest=72 # mean in lowest PRS decile
highest=86 # mean in lowest PRS decile
nSub=10000 # number of subjects to simulate

library(ggplot2)

subs=data.frame(matrix(0, ncol =4 , nrow =nSub))
colnames(subs)=c("Decile","zPRS","zWeight","Weight")
subs$zPRS=rnorm(nSub,sd=1)
subs$zWeight=subs$zPRS*sqrt(pVar)+rnorm(nSub,sd=sqrt((1-sqrt(pVar))))
t=sprintf("Weight in kilograms versus PRS, r = %.2f",r)
subs=subs[order(subs$zPRS),]
subs$Decile=10*seq(1,nSub)/nSub+.5

meanLow=mean(subs$zWeight[subs$Decile<1.5])
meanHigh=mean(subs$zWeight[subs$Decile>9.5])
subs$Weight=(subs$zWeight-meanLow)*(highest-lowest)/(meanHigh-meanLow)+lowest
WvPRS=ggplot(subs,aes(Decile,Weight)) + geom_point() + theme_bw() +
  scale_x_continuous(breaks = seq(1,10,by=1)) + 
  ggtitle(t)
WvPRS
cor(subs$Weight,subs$zPRS)
# just to check
cor(subs$Weight,subs$Decile)
mean(subs$Weight[subs$Decile<1.5])
mean(subs$Weight[subs$Decile>9.5])


