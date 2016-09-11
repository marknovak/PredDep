
# -------------------------------------------------------------------------
# Ensure parameter names match to permit vector of parameters to be passed,
# then fit models and print estimates
# -------------------------------------------------------------------------
Method = "BFGS"
REP = 100


print('*****************  Fitting density-independent model  ******************')
names(start_f) <- parnames(nll_f) <- paste0("f",1:PreySp)
fit_f = mle2(nll_f, start=start_f,  data=dat, method=Method, control=list(trace=1, REPORT=REP), hessian = SkipHess)
print(exp(coef(fit_f)))


print('*****************  Fitting Type II model  ******************')
names(start_t2) <- parnames(nll_t2) <- paste0("a",1:PreySp)
fit_t2 = mle2(nll_t2, start=start_t2, data=dat, method=Method, control=list(trace=1, REPORT=REP), hessian = SkipHess)
print(exp(coef(fit_t2)))


print('*****************  Fitting Ratio-depenent model  ******************')
names(start_rd) <- parnames(nll_rd) <- paste0("a",1:PreySp)
fit_rd = mle2(nll_rd, start=start_rd, data=dat, method=Method, control=list(trace=1, REPORT=REP), hessian = SkipHess)
print(exp(coef(fit_rd)))


print('*****************  Fitting Beddington-DeAngelis model ******************')
print('*****************   with intraspecific effects only   ******************')
names(start_bds) <- parnames(nll_bds) <- c(paste0("a",1:PreySp),'g1')
fit_bds = mle2(nll_bds, start=start_bds, data=dat, method=Method, control=list(trace=1, REPORT=REP, maxit=1000), hessian = SkipHess)
print(c(exp(coef(fit_bds)[1:PreySp]),coef(fit_bds)[-c(1:PreySp)]))


print('*****************  Fitting Beddington-DeAngelis model ******************')
print('*****************  with intra & interspecific effects ******************')
names(start_bdm) <- parnames(nll_bdm) <- c(paste0("a",1:PreySp),paste0("g",1:PredSp))
fit_bdm = mle2(nll_bdm, start=start_bdm, data=dat, method=Method, control=list(trace=1, REPORT=REP, maxit=1000), hessian = SkipHess)
print(c(exp(coef(fit_bdm)[1:PreySp]),coef(fit_bdm)[-c(1:PreySp)]))

