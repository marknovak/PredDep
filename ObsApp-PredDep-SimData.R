

# ##########################################################################
# Create multinomial data assuming Beddington-DeAngelis functional response
#    [Can be converted to Multispecies Type II by setting g's to zero]
# ##########################################################################
# --------------------------
# Specify measured variables
# --------------------------

# obs - total number of feeding observations on focal predator
# P -  densities of predators (across columns) in each surveyed site (across rows)
# N - densities of prey (across columns) in each surveyed site (across rows)
# h - handling times of prey (across columns) in each surveyed site (across rows)
obs = matrix(round(runif(1*Surv,obs.min,obs.max),0), nrow=1, ncol=Surv)
  rownames(obs) = 'obs'
P = matrix(round(runif(PredSp*Surv,P.min,P.max),0), nrow=PredSp, ncol=Surv)
  rownames(P) = paste0('P',1:PredSp)
N = matrix(round(runif(PreySp*Surv,N.min,N.max),0), nrow=PreySp, ncol=Surv)
  rownames(N) = paste0('N',1:PreySp)
h = matrix(round(runif(PreySp*Surv,h.min,h.max),0), nrow=PreySp, ncol=Surv)
  rownames(h) = paste0('h',1:PreySp)
  
# true.a - vector of focal predator's attack rates on each prey
# true.g - vector of focal predator's interference rates with each predator
true.a =rlnorm(PreySp*1,log(0.01/PreySp),1)
  names(true.a) = paste0('a',1:PreySp)
true.g = rlnorm(PredSp*1,log(0.1/PredSp),0.01)
  names(true.g) = paste0('g',1:PredSp)

# -------------
# Create dataset
# --------------
pi = t(t(true.a*N*h) / (1 + colSums(true.a*h*N) + colSums(true.g*P)))
p = rbind(p0=1-colSums(pi),pi)
  rownames(p) = paste0('p',0:PreySp)

n = t(rmultinomial(Surv,t(obs),t(p)))
  rownames(n) = paste0('n',0:PreySp)

# Associated 'true' feeding rates
true.f = t(t(true.a*N) / (1 + colSums(true.a*h*N) + colSums(true.g*P)))
  rownames(true.f) = paste0('f',1:PreySp)
  
# # Check to ensure no feeding rates exceed 1/h by random chance.
# anal_f = (dat$n[-1,]/dat$n[1,])*(1/(dat$h))
#   rownames(anal_f) = paste0('f',1:PreySp)
# any(true.f>1/h)
# any(anal_f>1/h)
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat = list(h=h, N=N, P=P, obs=obs, n=n)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

# ---------
# Plot data
# ---------
Nlim = c(0,max(N))
Rlim = c(0,max(t(N)/P[1,]))
Plim = c(0,max(P[1,]))
flim = c(0,max(true.f))
         
par(mfrow=c(1,3))
  matplot(t(N), t(true.f), xlim=Nlim, ylim=flim, pch=21, bg='grey', xlab=expression(N[i]), ylab=expression(f[i]))
  matplot(t(N)/P[1,], t(true.f), xlim=Rlim, ylim=flim, pch=21, bg='grey', xlab=expression(N[i]/P[1]), ylab=expression(f[i]))
  matplot(P[1,], t(true.f), xlim=Plim, ylim=flim, pch=21, bg='grey', xlab=expression(P[1]), ylab=expression(f[i]))
par(mfrow=c(1,1))
  # scatter3D(N,P,f,xlab='N',ylab='P',zlab='f',type='h',pch=19,theta=-45,phi=30,r=5)

################################################################################
################################################################################
################################################################################
