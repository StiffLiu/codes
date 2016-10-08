library("Ecdat")
?CPSch3
data(CPSch3)
dimnames(CPSch3)[[2]]

male.earnings = CPSch3[CPSch3[, 3]=="male", 2]
sqrt.male.earnings = sqrt(male.earnings)
log.male.earnings = log(male.earnings)

par(mfrow=c(2,2))
qqnorm(male.earnings, datax=T, main="untransformed")
qqnorm(sqrt.male.earnings, datax=T, main="square-root transformed")
qqnorm(log.male.earnings, datax=T, main="log-transformed")

par(mfrow=c(2,2))
boxplot(male.earnings, main="untransformed")
boxplot(sqrt.male.earnings, main="square-root transformed")
boxplot(log.male.earnings, main="log transformed")

par(mfrow=c(2,2))
plot(density(male.earnings), main="untransformed")
plot(density(sqrt.male.earnings), main="square-root transformed")
plot(density(log.male.earnings, main="log transformed")

library("MASS")
windows()
boxcox(male.earnings~1)
boxcox(male.earnings~1,lambda = seq(.3, .45, 1/100))
bc = boxcox(male.earnings~1,lambda = seq(.3, .45, by=1/100),interp=F)
ind = (bc$y==max(bc$y))
ind2 = (bc$y > max(bc$y) - qchisq(.95,df=1)/2)
bc$x[ind]
bc$x[ind2]

library("fGarch")
stdFit = sstdFit(male.earnings,hessian=T)
stdEst = stdFit$estimate
gedFit = sgedFit(male.earnings,hessian=T)
gedEst = gedFit$par

plot(density(male.earnings), xlim=c(-1,54), ylim=c(0.0, 0.06), main="")
par(new=T)
plot(density(rsstd(length(male.earnings), mean=stdEst[[1]], sd=stdEst[[2]], 
  nu=stdEst[[3]], xi=stdEst[[4]])), xlim=c(-1,54), ylim=c(0.0, 0.06), col="blue", main="")
par(new=T)
plot(density(rsged(length(male.earnings), mean=gedEst[[1]], sd=gedEst[[2]], 
  nu=gedEst[[3]], xi=gedEst[[4]])), xlim=c(-1,54), ylim=c(0.0, 0.06), col="red", main="")

data(Garch, package="Ecdat")
library("fGarch")
data(EuStockMarkets)
logDAX = diff(log(EuStockMarkets[, 1]))


loglikStd = function(x){
  f = -sum(log(dstd(logDAX, x[1], x[2], x[3])))
  f
}
start = c(mean(logDAX), sd(logDAX), 4)
fitStd = optim(start, loglikStd, method="L-BFGS-B", 
  lower=c(-.1, .001, 2.1),  upper=c(.1, 1, 20))
aicStd = 2 * fitStd$value + 2 * length(fitStd$par)
print(c("MLE=", round(fitStd$par, digits=5)))

loglikSstd = function(x){
  f = -sum(log(dsstd(logDAX, x[1], x[2], x[3], x[4])))
  f
}
start = c(mean(logDAX), sd(logDAX), 4, 1)
fitSstd = optim(start, loglikSstd, method="L-BFGS-B", 
  lower=c(-.1, .001, 2.1),  upper=c(.1, 1, 20))
aicSstd = 2 * fitSstd$value + 2 * length(fitSstd$par)
print(c("MLE=", round(fitSstd$par, digits=5)))

transLogDAX = qnorm(pstd(logDAX, mean=fitStd$par[1], sd=fitStd$par[2], nu=fitStd$par[3]))
plot(density(transLogDAX))
plot(density(logDAX))

library('fGarch')
gasFlow = read.csv('F:\\R\\SaDAfFE\\data\\GasFlowData.csv')
xlimit = c(3e5,13.5e5)
ylimit = c(0.0, 3.5e-6)
stdEst = sstdFit(gasFlow[[1]])$estimate
plot(density(gasFlow[[1]]), xlim=xlimit, ylim=ylimit, main="")
par(new=T)
plot(density(rsstd(length(gasFlow[[1]]), mean=stdEst[[1]], sd=stdEst[[2]], 
  nu=stdEst[[3]], xi=stdEst[[4]])), xlim=xlimit, ylim=ylimit, col="blue", main="")
  
start = c(1, 1)
loglikPos = function(theta) {
  - sum(log(dpois(y, lamda=theta[1] + theta[2] * x)))
}
mle = optim(start, loglikPos, hessian=T)
invFishInfo = solve(mle$hessian)
options(digits=4)
mle$par
mle$value
mle$convergence
sqrt(diag(invFishInfo))

library(evir)
library(fGarch)
data(bmw)
start_bmw = c(mean(bmw), sd(bmw), 4)
loglik_bmw = function(theta){
  -sum(log(dstd(bmw, mean=theta[1], sd=theta[2], nu=theta[3])))
}
mle_bmw = optim(start_bmw,  loglik_bmw, hessian=T)
fishInfo_bmw = solve(mle_bmw$hessian)


data(siemens)
n=length(siemens)
par(mfrow=c(3,2))
qqplot(siemens,qt(((1:n)-.5)/n,2),ylab="t(2) quantiles",
  xlab="data quantiles")
qqplot(siemens,qt(((1:n)-.5)/n,3),ylab="t(3) quantiles",
  xlab="data quantiles")
qqplot(siemens,qt(((1:n)-.5)/n,4),ylab="t(4) quantiles",
  xlab="data quantiles")
qqplot(siemens,qt(((1:n)-.5)/n,5),ylab="t(5) quantiles",
  xlab="data quantiles")
qqplot(siemens,qt(((1:n)-.5)/n,8),ylab="t(8) quantiles",
  xlab="data quantiles")
qqplot(siemens,qt(((1:n)-.5)/n,12),ylab="t(12) quantiles",
  xlab="data quantiles")

loglik_siemens = function(theta){
  -sum(log(dt(siemens, df=theta)))
}
mle_siemens = stdFit(siemens)