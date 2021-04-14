# Test DDbintau
# SRC 2021-02-16

rm(list = ls())
graphics.off()

library(stats)

source('DDbintau.r')

# save results
# Save post-processed DLM data
# X.raw is original series, X.z is z score, X.rawmean and X.rawsd are mean & sd for z score
# nl is DLM lags, delta is DLM delta, X.dlm is input to DLM, Tstep is time,
# local equilibria are X.eq and SD.eq, with ratio z.eq
# level (intercept) measures are level or stdlevel (for level/sd)
# Yyhat, B.ests, B.sd, errvar are DLM outputs
#save(nl,delta,X.raw,X.z,X.rawmean,X.rawsd,X.dlm,Tstep,X.eq,SD.eq,level,stdlevel,Z.eq,
#     Yyhat,B.ests,B.sd,errvar,file=Fname)
#
#load(file='Paul0811_nominalDLMresult.Rdata')
#fname=c('Paul0811_DDtauweightsv2.Rdata')
#load(file='Paul1315_nominalDLMresult.Rdata')
#fname=c('Paul1315_DDtauweightsv2.Rdata')
#
#load(file='Peter0811_nominalDLMresult.Rdata')
#fname=c('Peter0811_DDtauweightsv2.Rdata')
load(file='Peter1315_DLMresult.Rdata')
fname=c('DDbin9tau_Peter1315.Rdata')

#load(file='Tuesday1315_DLMresult.Rdata')
#fname=c('DDbin9tau_Tuesday1315.Rdata')

Xvar = na.omit(stdlevel) #[10:length(stdlevel)]
Tstep = Tstep[1:length(Xvar)]  # if Xvar = DLM output then trim Tstep

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(Tstep,Xvar,type='l',lwd=2,xlab='day',
     ylab='log10(phycocyanin)')
grid()

# Set up for binning method
# DDbins = function(Xvar,bw,ntau,nbin)
title = c('Pigment variate')
ntau = 9
nbin = 100
bw <- 0.3*sd(Xvar)  # try between 0.1*sd and 0.5*sd 
# run function
DDout = DDbins(Xvar,bw,ntau,nbin)
#outlist = list(D1s,D2s,sigmas,bin.mid,D1,D2,sigma)
# extract smoothed output
D1s = DDout[[1]]
D2s = DDout[[2]]
sigmas = DDout[[3]]
bin.mid= DDout[[4]]
# look at results of DDbintau
windows()
par(mfrow=c(3,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
plot(D1s$x,D1s$y,type='l',lwd=2,col='blue',xlab='state',ylab='D1')
abline(h=0,lty=2)
plot(sigmas$x,sigmas$y,type='l',lwd=2,col='red',xlab='state',ylab='sigma')
plot(D2s$x,D2s$y,type='l',lwd=2,col='red',xlab='state',ylab='D2')

# Find equilibria
sdrift = sign(D1s$y)
dsdrift = c(0,-diff(sdrift))
xeq = D1s$x[which(!dsdrift == 0)]
ixeq = which(!dsdrift == 0)  # indices of the equilibria

print('equilibria and their indices',quote=F)
print(xeq,quote=F)
print(ixeq,quote=F)

save(Tstep,Xvar,bin.mid,D1s,D2s,sigmas,xeq,file=fname)
