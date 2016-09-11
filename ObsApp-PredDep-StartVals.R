# ----------------------------------------------------------------
# Use analytical solutions to inform starting values
# Means correspond to MLE's yet often cause convergence failure.
# Therefore use offset multiplier of analytical estimates.
OM = 1
# ----------------------------------------------------------------

# Feeding rate
anal_f = t(t(dat$n[-1,])/colSums(dat$n))*(1/dat$h)
anal_f[is.na(anal_f)] = 0
if(PreySp==1){anal_f = t(anal_f)}
rownames(anal_f) = paste0('f',1:PreySp)
start_f = log(rowMeans(anal_f)*OM)


# Type II
anal_t2 = t(t(dat$n[-1,])/dat$n[1,])*(1/(dat$h*dat$N))
anal_t2[is.na(anal_t2)] = 0; anal_t2[is.infinite(anal_t2)] = 0;
if(PreySp==1){anal_t2 = t(anal_t2)}
rownames(anal_t2) = paste0('a',1:PreySp)
start_t2 = log(rowMeans(anal_t2)*OM)


# Ratio-dependent
anal_rd = t(t(dat$n[-1,])/dat$n[1,])*(dat$P[1,]/(dat$h*dat$N))
anal_rd[is.na(anal_rd)] = 0; anal_rd[is.infinite(anal_rd)] = 0;
if(PreySp==1){anal_rd = t(anal_rd)}
rownames(anal_rd) = paste0('a',1:PreySp)
start_rd = log(rowMeans(anal_rd)*OM)


# Beddington-DeAngelis - single and multi-predator versions
start_bds = c(log(apply(anal_t2,1,quantile,0.99,na.rm=TRUE)),'g1'=0.00)
start_bdm = c(log(apply(anal_t2,1,quantile,0.99,na.rm=TRUE)),rep(0.00,PredSp))
names(start_bdm)[-c(1:PreySp)] = paste0("g",1:PredSp)


# Hassell-Varley
start_hv = c(log(apply(anal_rd,1,quantile,0.99,na.rm=TRUE)),'m1'=log(0.5))


# #Type II w/ fixed exponent
# start_t3 = log(exp(start_t2)*2)

# # Crowley-Martin - single and multi-predator versions
# start_cms = c(log(apply(anal_t2,1,quantile,0.99,na.rm=TRUE)),'g1'=0.00)
# start_cmm = c(log(apply(anal_t2,1,quantile,0.99,na.rm=TRUE)),rep(0.00,PredSp))
# names(start_cmm)[-c(1:PreySp)] = paste0("g",1:PredSp)

