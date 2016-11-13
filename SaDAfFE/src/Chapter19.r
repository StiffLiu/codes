alpha = seq(0.025, 0.25, 0.005)
vaR = qnorm(1 - alpha, mean=-4000, sd=18000)
plot(alpha, vaR, type='l')

data(SP500, package="Ecdat")
n = length(SP500$r500)
years = 1981 + (1:n)* (1991.25-1981)/n
years = years[(n-999):n]
SPreturn = SP500$r500[(n-999):n]
S = 20000
alpha = 0.05
qAlpha = quantile(SPreturn, alpha)
vaR = -S * qAlpha
ES = -S * mean(SPreturn[SPreturn <= qAlpha])

library(MASS)
library(fGarch)
Tdist = fitdistr(SPreturn, densfun="t")
est = Tdist$estimate
n = length(SPreturn)
qqplot(SPreturn, qstd((1:n)/(n+1), mean=est['m'], sd=est['s'], nu=est['df']),
  main="t-probability plot", xlab="data", ylab="t-quantiles")
abline(lm(qstd(c(.25, .75), mean=est['m'], sd=est['s'], 
  nu=est['df']) ~ quantile(SPreturn, c(.25, .75))))

esT = function(mu, sl, df, alpha){
  lqt = qt(alpha, df=df)
  ldt = dt(lqt, df=df)
  as.numeric(-mu + sl * (ldt / alpha) * ((df + lqt^2) / (df - 1)))
}

esNorm = function(mu, sd, alpha){
  as.numeric(-mu + sd * dnorm(qnorm(alpha)) / alpha)
}

S * esT(est['m'], est['s'], est['df'], alpha)

bootStrapEst = function(dat, iterations=500, sampleSize=1000, alpha=0.05){
  npVaRs = rep(0, iterations)
  parVaRs = rep(0, iterations)
  npES = rep(0, iterations)
  parES = rep(0, iterations)

  for(i in 1:iterations){
    samples = sample(dat, sampleSize, replace=T)
    npVaRs[i] = quantile(samples, alpha)
    npES[i] = mean(samples[samples <= npVaRs[i]])
    fitTDist = fitdistr(samples, densfun="t")
    estT = fitTDist$estimate
    parVaRs[i] = estT['m'] + estT['s'] * qt(alpha, df=estT['df'])
    parES[i] = esT(estT['m'], estT['s'], estT['df'], alpha)
  }
  list(npVaRs=npVaRs, parVaRs=parVaRs, npES=npES, parES=parES)
}

attach(bootStrapEst(SPreturn))

# Non parametric booststrap estimations of value at risk
quantile(npVaRs, 0.95) * S
quantile(npVaRs, 0.05) * S
# Parametric booststrap estimations of value at risk
quantile(parVaRs, 0.95) * S
quantile(parVaRs, 0.05) * S
# Non parametric booststrap estimations of expected shortfall
quantile(npES, 0.95) * S
quantile(npES, 0.05) * S
# Parametric bootstrap estimations of expected shortfall
quantile(parES, 0.95) * S
quantile(parES, 0.05) * S

mu1 = mean(VaRs[1:50,1])
sd1 = sd(VaRs[1:50,1])
mu1 - qnorm(0.05) * sd1
mu1 + qnorm(0.05) * sd1

myFit = garchFit(formula=~garch(1, 1), data=SPreturn, cond.dist="std")
myFit
param = myFit@fit$par
pred = as.numeric(predict(myFit, n.ahead=1))
degreeFreedom = param['shape']
expectedShortFall = esT(mu=pred[1], sl=pred[3]/sqrt(degreeFreedom/(degreeFreedom - 2)),
  df=degreeFreedom, alpha) * S
valueAtRisk = qstd(alpha, mean=pred[1], sd=pred[3], nu=degreeFreedom) * S
plot(years, myFit@sigma.t, type='l')
# abline(h=sd(SPreturn))
abline(h = Tdist$estimate[2]*sqrt((Tdist$estimate[3])/(Tdist$estimate[3] - 2)), lty=5)
points(years[n] + 1/365, pred[3], cex=4, pch='*')
legend("topright", c("conditional SD", "marginal SD", "next day's conditional SD"),
  lty=c(1, 5, NA), pch=c(NA, NA, '*'), pt.cex=2.5)

library(MASS)
data(CRSPday, package='Ecdat')
colNames = c('ge', 'ibm', 'mobil')
dat = CRSPday[, colNames]

library('mnormt')
dfs = seq(5.7, 5.726, length=500)
logLiks = rep(0, length(dfs))
for(i in 1:length(dfs)){
  ret = cov.trob(dat, nu=dfs[i])
  logLiks[i] = sum(dmt(dat, mean=ret$center, ret$cov, dfs[i], log=T))
}
plot(dfs, logLiks)
df = dfs[(1:length(logLiks))[logLiks == max(logLiks)]]
portfolio = cov.trob(dat, nu=df)

wghts = c(1/3, 1/3, 1/3)
mu = sum(portfolio$center * wghts)
variance = t(wghts) %*% portfolio$cov %*% wghts
sl = sqrt((df - 2)/df * variance)
valueAtRisk = -S * (qt(alpha, df=df) * sl + mu)
expectedShortFall = esT(mu, sl, df, alpha) * S

SSPR = sort(SPreturn)
plotPolyTail = function(m){
  selected = 1:m
  logSSPR = log(-SSPR[selected])
  logIndex = log(selected / length(SSPR))
  fit = coef(lm(logSSPR~logIndex))
  slope = fit['logIndex']
  intercept = fit['(Intercept)']
  plot(logIndex, logSSPR, xlab="log(k/n)", ylab=expression(log(-R[(k)])),
    main=paste("m=", m, ", slope=", slope, ", a=", -1/slope))
  abline(coef=fit, lwd=2)
  -1/slope
}

par(mfrow=c(2,2))
plotPolyTail(50)
polyTailRegIndex=100
regressEst = plotPolyTail(polyTailRegIndex)
plotPolyTail(200)
plotPolyTail(300)

nc = 20:250
estA = rep(0, length(nc))
# SSPR is sorted.
for(i in 1:length(nc)){
  estA[i] = nc[i]/sum(log(SSPR[1:nc[i]]/SSPR[nc[i]]))
}
par(mfrow=c(1, 3))
plot(nc, estA, ylab="Hill estimator", main="(a)")
plot(nc, estA, ylab="Hill estimator", main="(b)", xlim=c(25, 120))
plot(nc, estA, ylab="Hill estimator", main="(c)", xlim=c(60, 100))

# from the above plot, we take 2.2 for Hill Estimator
hillEstimateIndex = (1:length(estA))[abs(estA - 2.2) == min(abs(estA - 2.2))]
hillEstimate = estA[hillEstimateIndex]

library(MASS)
Tdist = fitdistr(SPreturn, densfun="t")$estimate
alphas = seq(0.001, 0.25, length=1000)
varNorm = -S *(mean(SPreturn) + qnorm(alphas) * sd(SPreturn))
varT = -S * (qt(alphas, df=Tdist['df']) * Tdist['s'] + Tdist['m'])
alpha = 0.1
vaR = quantile(SSPR, alpha)
varPolyTailReg = -S * vaR * ((alpha/alphas)^(1/regressEst))
varPolyTailHill = -S * vaR * ((alpha/alphas)^(1/hillEstimate))

plot(alphas, varPolyTailReg, xlab=expression(alpha), 
  ylab=expression(VaR(alpha)), log='x', type='l', lty=1)
lines(alphas, varPolyTailHill, log='x', type='l', lty=2)
lines(alphas, varT, log='x', type='l', lty=3)
lines(alphas, varNorm, log='x', type='l', lty=4)
legend("topright", c=("polynomial tail: regression", "polynomial tail: Hill",
  "t", "normal"), lty=c(1, 2, 3, 4))

library(fEcofin)
library(mnormt)
Berndt = berndtInvest[,5:6]
logLikhood = function(theta){
  A = matrix(data=c(theta[3], theta[4], 0, theta[5]), nrow=2, byrow=T)
  -sum(log(dmt(Berndt, mean=c(theta[1], theta[2]), S=t(A) %*% A, df=theta[6])))
}

estimate = function(dat, alpha=0.05, wgts=c(0.3, 0.7)){
  initA0 = chol(cov(dat))
  initA = initA0[initA0 != 0]
  initMeans = apply(dat, 2, mean)
  fitMvt = optim(c(initMeans, initA, 5), logLikhood)
  means = fitMvt$par[1:length(initMeans)]
  covA = matrix(rep(0, nrow(initA0) * ncol(initA0)), nrow=nrow(initA0), ncol=ncol(initA0))
  covA[!lower.tri(covA)] = c(fitMvt$par[(length(initMeans) + 1):(length(initMeans) + length(initA))])
  degreeFreedom = fitMvt$par[length(initA) + length(initMeans) + 1]
  covariance = t(covA) %*% covA
  portMean = sum(wgts * means)
  portVar = wgts %*% covariance %*% wgts
  portScale = sqrt(portVar * (degreeFreedom - 2) / degreeFreedom)
  expectedShortFall = esT(portMean, portScale, degreeFreedom, alpha)
  valueAtRisk = qt(alpha, df=degreeFreedom) * portScale + portMean
  c(-valueAtRisk, expectedShortFall, degreeFreedom)
}

estimate(Berndt)

iterations = 250
estimations = matrix(nrow=iterations, ncol=3)
colnames(estimations) = c('valueAtRisk', 'expectedShortFall', 'degreeFreedom')
sampleSize = 60

for(i in 1:iterations){
  samples = Berndt[sample(1:nrow(Berndt), sampleSize, replace=T), ]
  estimations[i, ] = estimate(samples)
}

SDEC = sort(Berndt[, 'DEC'])
nc = 10:48
estA = rep(0, length(nc))
# SEC is sorted.
for(i in 1:length(nc)){
  estA[i] = nc[i]/sum(log(SDEC[1:nc[i]]/SDEC[nc[i]]))
}
plot(nc, estA, ylab="Hill estimator")

library(fEcofin)
library(MASS)
alpha = 0.01
bmwDat = bmwRet[, 2]
vaRnp = quantile(bmwDat, alpha)
esnp = mean(bmwDat[bmwDat <= vaRnp])

vaRNorm = qnorm(alpha, mean=mean(bmwDat), sd=sd(bmwDat))
esNorm = esNorm(mean(bmwDat), sd(bmwDat), alpha)

fitTDist = fitdistr(bmwDat, densfun="t")
estT = fitTDist$estimate
vaRT = estT['m'] + estT['s'] * qt(alpha, df=estT['df'])
expSfT = esT(estT['m'], estT['s'], estT['df'], alpha)

closingPrice = msft.dat[, 'Close']
dailyReturn = diff(log(closingPrice))
attach(bootStrapEst(dailyReturn, sampleSize=150))

dat = read.csv("F:/R/SaDAfFE/data/Stock_FX_Bond.csv",header=T)
prices = as.matrix(dat[1:500,c(3,5,7,9,11)])
returns = apply(log(prices), 2, diff)
rownames(returns) = NULL
mus = apply(returns, 2, mean)
cov(returns)
50 / length(prices[500, ]) / prices[500, ]
weights = (1/prices[500, ])/sum(1/prices[500, ])
portReturns = diff(log(prices %*% weights))
# Assuming normal distrutions
fitNorm = fitdistr(portReturns, "normal")$estimate
alpha = 0.01
# Value At Risk:
-qnorm(alpha, mean=fitNorm["mean"], sd=fitNorm["sd"])
# Expected shortfall:
esNorm(fitNorm["mean"], fitNorm["sd"], alpha)
# Assuming daily returns are uncorrelated.
# Then five day return are N(5*mean, sd*sqrt(5))
# Value At Risk:
-qnorm(alpha, mean=5*fitNorm["mean"], sd=fitNorm["sd"]*sqrt(5))
# Expected shortfall:
esNorm(5*fitNorm["mean"], fitNorm["sd"]*sqrt(5), alpha)