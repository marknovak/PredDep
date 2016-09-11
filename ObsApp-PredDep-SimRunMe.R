###############################################################################
###############################################################################
###############################################################################
rm(list=ls())

set.seed(100)

library(bbmle)
library(mc2d) # dmultinomial (vectorized version of dmultinom)

# -------------------
# Specify system info
# -------------------
Surv = 100
PreySp = 3
PredSp = 2
obs.min = 500
obs.max = 1000
P.min = 10
P.max = 100
N.min = 10
N.max = 200
h.min = 2
h.max = 5

Abund.scale = 1

# ------------
# Specify data
# ------------
source('ObsApp-PredDep-SimData.R')

# Rescale abundances to avoid numerical issues
dat$N = dat$N/Abund.scale
dat$P = dat$P/Abund.scale

# ~~~~~~~~~~~~~~~~~
rm(N, P, h, n, obs)
# ~~~~~~~~~~~~~~~~~
# -------------------
# Load and fit models
# -------------------
source('ObsApp-PredDep-Models.R')
source('ObsApp-PredDep-StartVals.R')


# -------------------------------------------------------
# Plot analytical estimates, truths, and starting values
# -------------------------------------------------------
# Feeding rate
par(pty='s')
lims=c(0,max(c(true.f,anal_f)))
matplot(t(true.f), t(anal_f), pch=21, bg='grey', 
        ylab='Est. feeding rate', xlab='True feeding rate', xlim=lims, ylim=lims)
abline(0,1)

# Type II
par(mfrow=c(1,PreySp),pty='m')
for(i in 1:PreySp){
  hist(anal_t2[i,], ylab='', 
       xlab='Est. Type II\nattack rate', col='grey', main='', breaks=20)
  abline(v=true.a[i],lwd=4)
  abline(v=start_t2[i],lwd=4,lty=2)
  legend('topright',bty='n',inset=0.1,c('Start','True'),lwd=2,lty=c(2,1))}

# Ratio-dependent
par(mfrow=c(1,PreySp),pty='m')
for(i in 1:PreySp){
  hist(anal_rd[i,], ylab='', 
       xlab='Est. ratio-dependent\nattack rate', col='grey', main='', breaks=20)
  abline(v=start_rd[i],lwd=4,lty=2)}
par(mfrow=c(1,1))


SkipHess=TRUE  # Skip Hessian if interested in point estimates only
source('ObsApp-PredDep-FitModels.R')       

# --------------
# Compare models
# --------------
# Notes: 
# - For the simulated dataset, density-independent model will only do best with (a) no variation in (prey) abundances *and* no variation in handling times across surveys.  If either is the case, then the Type II will be the best-performing when no predator interference occurs.

models=c('Density-independent','Type II','Ratio-dependent','Beddington-DeAngelis - intra', 'Beddington-DeAngelis - inter')

tab_AIC = AICtab(fit_f,fit_t2,fit_rd,fit_bds,fit_bdm, base=TRUE, logLik=FALSE, weights=TRUE, mnames=models)
tab_AICc = AICctab(fit_f,fit_t2,fit_rd,fit_bds,fit_bdm, base=TRUE, logLik=FALSE, weights=TRUE, mnames=models, nobs=ncol(dat$n))

# ----------------------------------------------
# Compare estimates to truth (for simulated data)
# ----------------------------------------------
tab_ests = array(NA,dim=c(length(models),PreySp+PredSp))
  colnames(tab_ests) = c(names(true.a),names(true.g))
  rownames(tab_ests) = c('True value', models[-1])
  tab_ests[1,] = c(true.a, true.g)
  tab_ests[2,1:length(coef(fit_t2))] = exp(coef(fit_t2)) / Abund.scale
  tab_ests[3,1:length(coef(fit_rd))] = exp(coef(fit_rd)) / Abund.scale
  tab_ests[4,1:length(coef(fit_bds))] = exp(coef(fit_bds)) / Abund.scale
  tab_ests[5,1:length(coef(fit_bdm))] = exp(coef(fit_bdm)) / Abund.scale
tab_ests = format(tab_ests,digits=3)
tab_ests[grep('NA',tab_ests)] = '-'
tab_ests = data.frame(tab_ests)


tab_AIC = data.frame(print(tab_AIC))
tab_AICc = data.frame(print(tab_AICc))
tab_ests


par(mfrow=c(1,2),pty='s')
 plot(true.a,exp(coef(fit_t2)),type='n',xlab='True',ylab='Est. Type II')
  abline(0,1)
  points(true.a,exp(coef(fit_t2)),pch=21,bg='grey')

plot(true.a,exp(coef(fit_bdm)[1:PreySp]),type='n',xlab='True',ylab='Est. BD')
  abline(0,1)
  points(true.a,exp(coef(fit_bdm))[1:PreySp],pch=21,bg='grey')
par(pty='m')


################################################################################
################################################################################
################################################################################