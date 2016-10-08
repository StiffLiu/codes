library('fEcofin')
library('fUtilities')
library('fGarch')

n = dim(bmwRet)[1]
kurt = kurtosis(bmwRet[,2], method="moment")
skew = skewness(bmwRet[,2], method="moment")
fit_skewt = sstdFit(bmwRet[,2])
q.grid = (1:n) / (n + 1)
stdEst = fit_skewt$estimate
qqplot(bmwRet[,2], qsstd(q.grid, stdEst[1], stdEst[2],
  stdEst[3], stdEst[4]), ylab="skewed-t quantiles")
skewness(rsstd(10000, stdEst[1], stdEst[2],
  stdEst[3], stdEst[4]), method="moment")
  
quantKurt.function = function(f, param, p1=0.025, p2=0.25){
  if (p1 >= 0.5) stop("p1 should be less than 0.5")
  if (p2 >= 0.5) stop("p2 should be less than 0.5")
  if (p1 >= p2) stop("p1 should be greater than p2")
  Q = f(c(p1, p2, 1 - p2, 1 - p1), param)
  k = (Q[4] - Q[1]) / (Q[3] - Q[2])
  k
}
quantKurt.numeric = function(y, p1=0.025, p2=0.25){
  if (p1 >= 0.5) stop("p1 should be less than 0.5")
  if (p2 >= 0.5) stop("p2 should be less than 0.5")
  if (p1 >= p2) stop("p1 should be greater than p2")
  Q = quantile(y, c(p1, p2, 1 - p2, 1 - p1))
  k = (Q[4] - Q[1]) / (Q[3] - Q[2])
  k
}
tQuantKurt = function(dfs){
  len = length(dfs)
  kurts = rep(0, len)
  for (i in 1:len){
    kurts[i] = quantKurt.function(qt, dfs[i])
  }
  kurts
}
dfs = seq(1, 10, by=0.02)
kurts = tQuantKurt(dfs)
plot(dfs, kurts)

nboot = 5000
n = 100
ModelFree_kurt = rep(0, nboot)
ModelBased_kurt = rep(0, nboot)

set.seed("5640")
for(i in 1:nboot){
  samp_ModelFree = sample(bmwRet[,2], n, replace = TRUE)
  samp_ModelBased = rsstd(n, stdEst[1], stdEst[2], stdEst[3], stdEst[4])
  ModelFree_kurt[i] = quantKurt.numeric(samp_ModelFree)
  ModelBased_kurt[i] = quantKurt.numeric(samp_ModelBased)
}
plot(density(ModelFree_kurt), xlim=c(1.9, 10.2), ylim=c(0, 0.58))
par(new=TRUE)
plot(density(ModelBased_kurt), xlim=c(1.9, 10.2), ylim=c(0, 0.58), col='red')
alpha = 0.90
quantile(ModelFree_kurt, c(0.5 - alpha / 2, 0.5 + alpha / 2))
quantile(ModelBased_kurt, c(0.5 - alpha / 2, 0.5 + alpha / 2))

