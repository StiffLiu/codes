library('fEcofin')
Berndt = as.matrix(berndtInvest[, 2:5])
cov(Berndt)
cor(Berndt)
pairs(Berndt)

library('MASS')
library('mnormt')
df = seq(2.5, 8, 0.01)
n = length(df)
loglik_max = rep(0, n)
for(i in 1:n) {
  fit = cov.trob(Berndt, nu=df[i])
  mu = as.vector(fit$center)
  sigma = matrix(fit$cov, nrow=4)
  loglik_max[i] = sum(log(dmt(Berndt, mean=fit$center, S=fit$cov, df=df[i])))
}



dfMax = df[loglik_max == max(loglik_max)]
dfConfInterval = df[loglik_max > max(loglik_max) - qchisq(0.95, df=1)/2]

library("MASS")
library("fGarch")
par(mfrow=c(1,4))
N = 2500
nu = 3
numSeed = 5640

means = c(0.001, 0.002)
sigma = matrix(c(0.1, 0.03, 0.03, 0.15), nrow=2)
weight = matrix(c(0.5, 0.5), nrow=1)
variance = weight %*% sigma %*% t(weight)
average = weight %*% matrix(means, nrow=2)
tVariables1 = rstd(N, mean=average, sd=variance, nu=5)
set.seed(numSeed)
chisqVars = sqrt(5 / rchisq(N, df=5))
tVariablesTemp = mvrnorm(N, mu=means, Sigma=sigma) * cbind(chisqVars, chisqVars)
tVariables2 = tVariables2[,1] + tVariables2[,2]
qqplot(tVariables1, tVariables2)

set.seed(numSeed)
cov = matrix(c(1, .8, .8, 1), nrow=2)
x = mvrnorm(N, mu=c(0, 0), Sigma=cov)
w = sqrt(nu / rchisq(N, df=nu))
x = x * cbind(w, w)
plot(x, main="(a)")

set.seed(numSeed)
cov = matrix(c(1, .8, .8, 1), nrow=2)
x = mvrnorm(N, mu=c(0, 0), Sigma=cov)
w1 = sqrt(nu/rchisq(N, df=nu))
w2 = sqrt(nu/rchisq(N, df=nu))
x = x * cbind(w1, w2)
plot(x, main="(c)")

set.seed(numSeed)
cov = matrix(c(1, 0, 0, 1), nrow=2)
x = mvrnorm(N, mu=c(0, 0), Sigma=cov)
w1 = sqrt(nu/rchisq(N, df=nu))
w2 = sqrt(nu/rchisq(N, df=nu))
x  = x * cbind(w1, w2)
plot(x, main="(c)")

set.seed(numSeed)
cov = matrix(c(1, 0, 0, 1), nrow=2)
x = mvrnorm(N, mu=c(0, 0), Sigma=cov)
w = sqrt(nu/rchisq(N, df=nu))
x = x * cbind(w, w)
plot(x, main="(d)")

library('mnormt')
data(CRSPday, package='Ecdat')
Y = CRSPday[, c(5,7)]
loglik = function(par) {
  mu = par[1:2]
  A = matrix(c(par[3], par[4], 0, par[5]), nrow=2, byrow=T)
  scale_matrix = t(A) %*% A
  df = par[6]
  f = -sum(log(dmt(Y, mean=mu, S=scale_matrix, df=df)))
  f
}
A = chol(cov(Y))
start = as.vector(c(apply(Y, 2, mean), A[1, 1], A[1, 2], A[2,2], 4))
fit_mvt = optim(start, loglik, method="L-BFGS-B", lower=c(-.02, -.02, -.1, -.1, -.1, 2), upper=c(.02, .02, .1, .1, .1, 15), hessian=T)
