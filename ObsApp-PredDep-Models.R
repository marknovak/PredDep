
# -------------------------------
# Define negative log likelihoods
# -------------------------------
# Feeding rate
nll_f = function(parms){
  f = exp(parms)
  if (any(0 > f*h | f*h > 1)){return(NA)}
  pi = f*h
  if(PreySp==1){  p = rbind(p0=1-sum(pi),pi) }
  if(PreySp>1){   p = rbind(p0=1-colSums(pi),pi) }
  if (any(0 > p | p > 1)){return(NA) }
  -sum(dmultinomial(t(n), prob=t(p), log=TRUE))
}

# Type II
nll_t2 = function(parms){
  a = exp(parms)
  # print(paste("a =",a))
  if(any(a > 999)){return(NA)}
  if(PreySp==1){
    pi = a*N*h / (1 + a*h*N)
    p = rbind(p0=1-pi,pi)}
  if(PreySp>1){
    pi = t(t(a*N*h) / (1 + colSums(a*h*N)))
    p = rbind(p0=1-colSums(pi),pi) }
  if (any(0 > p | p > 1)){return(NA) }
  -sum(dmultinomial(t(n), prob=t(p), log=TRUE))
}

# Ratio-dependent
nll_rd = function(parms){
  a = exp(parms)
  # print(paste("a =",a))
  if(any(a > 999)){return(NA)}
  if(PreySp==1){
    pi = a*N*h / (P[1,] + a*h*N)
    p = rbind(p0=1-pi,pi) }
  if(PreySp>1){
    pi = t(t(a*N*h) / (P[1,] + colSums(a*h*N)))
    p = rbind(p0=1-colSums(pi),pi) }
  if (any(0 > p | p > 1)){return(NA) }
  -sum(dmultinomial(t(n), prob=t(p), log=TRUE))
}

# Beddington-DeAngelis - intraspecific predator effects only
nll_bds = function(parms){
  a = exp(parms[1:PreySp])
  g = (parms[-c(1:PreySp)])
  # print(c(a,g))
  if(any(c(a,g) > 999)){return(NA)}
  if(PreySp==1){
    pi = a*N*h / (1 + a*h*N + g*P[1,]) # only difference to nll_bdm
    p = rbind(p0=1-pi,pi)}
  if(PreySp>1){
    pi = t(t(a*N*h) / (1 + colSums(a*h*N) + g*P[1,])) # only difference to nll_bdm
    p = rbind(p0=1-colSums(pi),pi)}
  # print(p)
  if (any(0 > p | p > 1)){return(NA) }
  -sum(dmultinomial(t(n), prob=t(p), log=TRUE))
}

# Beddington-DeAngelis - intra + interspecific predator effects
nll_bdm = function(parms){
  a = exp(parms[1:PreySp])
  g = (parms[-c(1:PreySp)])
  # print(exp(parms))
  if(any(c(a,g) > 999)){return(NA)}
  if(PreySp==1){
    pi = a*N*h / (1 + a*h*N + colSums(g*P))
    p = rbind(p0=1-pi,pi)}
  if(PreySp>1){
    pi = t(t(a*N*h) / (1 + colSums(a*h*N) + colSums(g*P))) 
    p = rbind(p0=1-colSums(pi),pi) }
  if (any(0 > p | p > 1)){return(NA) }
  -sum(dmultinomial(t(n), prob=t(p), log=TRUE))
}

# Hassell-Varley
nll_hv = function(parms){
  a = exp(parms[1:PreySp])
  m = exp(parms[-c(1:PreySp)])
  # print(paste(c(a,m)))
  if(any(c(a,m) == 0)){return(NA)}
  if(any(abs(c(a,m)) > 999)){return(NA)}
  if(PreySp==1){
    pi = a*N*h / (P[1,]^m + a*h*N)
    p = rbind(p0=1-pi,pi) }
  if(PreySp>1){
    pi = t(t(a*N*h) / (P[1,]^m + colSums(a*h*N)))
    p = rbind(p0=1-colSums(pi),pi) }
  # print(p)
  if (any(0 > p | p > 1)){return(NA) }
  -sum(dmultinomial(t(n), prob=t(p), log=TRUE))
}
