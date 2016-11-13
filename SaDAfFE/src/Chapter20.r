library(fGarch) #for function dstd

iterations = 30000
results = matrix(ncol=3, nrow=iterations)
numSamples = 500

set.seed(19485)
y = rstd(numSamples,mean=0, sd=1, nu=5)
results[1, ] = c(runif(1,.1,5), runif(1,.1,5), runif(1,2.5,15))

mus = seq(-3, 3, 0.01)
ps1 = rep(0, length(mus))
for(i in 1:length(ps1)){
  ps1[i] = sum(log(dstd(y, mean=mus[i], sd=1, nu=5)))
}
plot(mus, ps1)
mus[ps1 == max(ps1)]

sds = seq(0.01, 3, 0.01)
ps2 = rep(0, length(sds))
for(i in 1:length(ps2)){
  ps2[i] = sum(log(dstd(y, mean=0, sd=sds[i], nu=5)))
}
plot(sds, ps2)
sds[ps2 == max(ps2)]

dfs = seq(2.01, 8, 0.01)
ps3 = rep(0, length(dfs))
for(i in 1:length(ps3)){
  ps3[i] = sum(log(dstd(y, mean=0, sd=1, nu=dfs[i])))
}
plot(dfs, ps3)
dfs[ps3 == max(ps3)]

sum(log(dstd(y, mean=mus[ps1 == max(ps1)], sd=sds[ps2 == max(ps2)], nu=dfs[ps3 == max(ps3)])))

# A simplified version of Metropolis-Hastings algorithm?
# Assuming uniform priors
for(i in 2:iterations){

  # assuming guass proposal density based on previous assumption
  newMu = rnorm(1, mean=results[i - 1, 1])
  newScale = rnorm(1, mean=results[i - 1, 2])
  newDf = rnorm(1, mean=results[i - 1, 3])

  # assuming guass proposal density
  # newMu = rnorm(1, mean=0)
  # newScale = rnorm(1, mean=1)
  # newDf = rnorm(1, mean=5)
  
  # assuming uniform proposal density
  #newMu = runif(1, min=-0.5, max=0.5)
  #newScale = runif(1, min=0.0, max=2)
  #newDf = runif(1, min=2.5, max=7.5)
  
  prob = sum(log(dstd(y, results[i-1,1], results[i-1,2], results[i-1,3])))
  logr = sum(log(dstd(y, newMu, results[i-1,2], results[i-1,3]))) - prob
  # If "logr" is not a value, that means we have an almost zero probability for one value in the sample "y"
  #   so we need to use previous value.
  #   Note, this is very important, previously without this condition the estimation diverges very far
  #     from the true value and it's taking me a long time to find out the reason.
  #   
  # Without the:
  #   runif(1) > exp(logr)
  # the estimation will always go to the region with higher probability,
  #   which will leave small variations of the estimations, as we're estimations from
  #   the samples of true distribution, less variation means more bias.
  if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newMu = results[i-1,1]
  logr = sum(log(dstd(y, newMu, newScale, results[i-1,3]))) - prob
  if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newScale = results[i-1,2]
  logr = sum(log(dstd(y, newMu, newScale, newDf))) - prob
  if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newDf = results[i-1,3]
  results[i,] = c(newMu, newScale, newDf)
  
  # The theory behind the above code is:
  #   Let
  #     pi`: be our prior distribution for the parameters: mu, sd, df, 
  #          in our above code we're assumming uniform distribution
  #     pi: posterior distribution for the parameters: mu, sd, df
  #     pr: proposal density, in our code, it's gauss or uniform distribution, which is symmetric
  #          that is pr(a|b) == pr(b|a)
  #     r = (pi(theta`|Y) * pr(theta`|theta) * pi`(theta`)) / (pi(theta|Y) * pr(theta|theta`) * pi`(theta))
  #     theta`: estimation of the parameters: mu, sd and df.
  #     theta: estimation of the parameters: mu, sd and df, in the previous iteration
  #   if "r" is greater than or equal to 1, choose theta` for this iteration.
  #   if "r" is less than 1, choose theta` with probability "r"
  #   
  #   From the above we have pi`(theta`) = pi`(theta)(this equation does not required anyway, pr(theta`|theta) = pr(theta|theta`)
  #      so r = pi(theta`|Y) / pi(theta|Y)
  #   And again, by Bayese theory:
  #        pi(theta|Y) = pi`(theta) * pi(Y|theta) / intg(pi`(theta) * pi(Y|theta)), 
  #               where intg is the integration of the function pi`(theta) * pi(Y|theta), and is a constant
  #   Similiarly:
  #        pi(theta`|Y) = pi`(theta`) * pi(Y|theta`) / intg(pi`(theta`) * pi(Y|theta`)), 
  #               where intg is the integration of the function pi`(theta`) * pi(Y|theta`), and is a constant
  #   And we have:
  #       r = pi(Y|theta`) / pi(Y|theta)
  #   Then: logr = log(pi(Y|theta`)) - log(pi(Y|theta))
}

# Hierarchical Priors
data('midcapD.ts', package='fEcofin')
x  = 100*as.matrix(midcapD.ts[,-c(1,22)])
m  = ncol(x)
n  = nrow(x)
k  = 100
x1 = x[1:k, ]
x2 = x[(k+1):n, ]
iterations = 5000
means = matrix(nrow=iterations, ncol=m)
means[1, ] = apply(x1, 2, mean)
alphas = rep(0.001, iterations)
eps = rep(sd(means[1, ]), iterations)
sds = rep(mean(apply(x1, 2, sd)), iterations)

#1 derive alphas based on means
for(i in 2:iterations){
  for(j in 1:m){
    # This loop does not take into consideration the assumption 
	# that the means are normally distribution
    oldMu = means[i - 1, j]
    oldEps = eps[if(j == 1) i - 1 else i]
    newMu = runif(1, min=oldMu/2, max=3*oldMu/2)
    newEps = rnorm(1, mean=oldEps)	
    prob = sum(log(dnorm(x1[, j], mean=oldMu, sd=oldEps)))
	
    logr = sum(log(dnorm(x1[, j], mean=newMu, sd=oldEps))) - prob
    if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newMu = oldMu

    logr = sum(log(dnorm(x1[, j], mean=newMu, sd=newEps))) - prob
    if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newEps = oldEps

    means[i, j] = newMu
    eps[i] = newEps
  }
  
  oldAlpha = alphas[i - 1]
  oldSd = sds[i - 1]
  newAlpha = rnorm(1, mean=oldAlpha/2, 3*oldAlpha/2)
  newSd = rnorm(1, mean=oldSd)
  prob = sum(log(dnorm(means[i, ], mean=oldAlpha, sd=oldSd)))

  logr = sum(log(dnorm(means[i, ], mean=newAlpha, sd=oldSd))) - prob
  if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newAlpha = oldAlpha

  logr = sum(log(dnorm(means[i, ], mean=newAlpha, sd=newSd))) - prob
  if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newSd = oldSd

  alphas[i] = newAlpha
  sds[i] = newSd
}

#2 derive means based on alpha
for(i in 2:iterations){
  oldAlpha = alphas[i - 1]
  oldSd = sds[i - 1]
  newAlpha = rnorm(1, mean=oldAlpha)
  newSd = rnorm(1, mean=oldSd)
  prob = sum(log(dnorm(means[i - 1, ], mean=oldAlpha, sd=oldSd)))

  logr = sum(log(dnorm(means[i - 1, ], mean=newAlpha, sd=oldSd))) - prob
  if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newAlpha = oldAlpha

  logr = sum(log(dnorm(means[i - 1, ], mean=newAlpha, sd=newSd))) - prob
  if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newSd = oldSd

  alphas[i] = newAlpha
  sds[i] = newSd

  for(j in 1:m){
    oldMu = means[i - 1, j]
    oldEps = eps[if(j == 1) i - 1 else i]
    newMu = rnorm(1, mean=newAlpha, sd=newSd)
    newEps = rnorm(1, mean=oldEps)	
    prob = sum(log(dnorm(x1[, j], mean=oldMu, sd=oldEps)))
	
    logr = sum(log(dnorm(x1[, j], mean=newMu, sd=oldEps))) - prob
    if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newMu = oldMu

    logr = sum(log(dnorm(x1[, j], mean=newMu, sd=newEps))) - prob
    if (is.nan(logr) || (logr < 0 && runif(1) > exp(logr))) newEps = oldEps

    means[i, j] = newMu
    eps[i] = newEps
  }
}

mu1 = apply(x1,2,mean)
mu2 = apply(x2,2,mean)
meanOfMeans = apply(means, 2, mean)
par(mfrow=c(1,2))
plot(c(rep(1,m),rep(2,m)),c(mu1,mu2),
   xlab="estimate                         target",ylab="mean",
   main="sample means",
   ylim=c(-.4,.8),axes=F)
axis(2)
axis(1,labels=F,tick=T,lwd.tick=0)
for (i in 1:m){ lines(1:2,c(mu1[i],mu2[i]),col=1) }

plot(c(rep(1,m),rep(2,m)),c(meanOfMeans,mu2),
   xlab="estimate                         target",ylab="mean",
   main="Bayes",
   ylim=c(-.4,.8),axes=F)
axis(2)
axis(1,labels=F,tick=T,lwd.tick=0)
for (i in 1:m){ lines(1:2,c(meanOfMeans[i],mu2[i]),col=1) }

library(R2WinBUGS)
data(CRSPmon, package="Ecdat")
ibm = CRSPmon[,2]
y = as.numeric(ibm)
N = length(y)
ibm_data = list("y", "N")

tbrate_t_bugs_code="
model{
  for(i in 1:N){
    y[i] ~ dt(mu, tau, nu)
  }
  
  # Prior distribution
  mu ~ dnorm(0.0, 1.0E-6)
  tau ~ dgamma(0.1, 0.01)
  sigma <- 1 / sqrt(tau)
  nu ~ dunif(2, 50)
}
"
tbrate_t_bugs_filename = "Tbrate_t.bug"
tbrate_t_bugs_file = file(tbrate_t_bugs_filename)
writeLines(tbrate_t_bugs_code, tbrate_t_bugs_file)
close(tbrate_t_bugs_file)

inits = function(){list(mu=rnorm(1, 0, .3),
  tau=runif(1, 1, 10), nu=runif(1, 1, 30))}
univt.mcmc = bugs(ibm_data, inits, model.file=tbrate_t_bugs_filename,
  parameters=c("mu", "tau", "nu", "sigma"),
  n.chains=3, n.iter=2600, n.burnin=100, n.thin=1,
  bugs.directory="E:/windows/programs/WinBUGS14",
  codaPkg=F, bugs.seed=5640)
file.remove(tbrate_t_bugs_filename)
print(univt.mcmc, digits=4)
# plot(univt.mcmc)

chains.mu = matrix(univt.mcmc$sims.array[,,'mu'], ncol=1)
chains.tau = matrix(univt.mcmc$sims.array[,,'tau'], ncol=1)
chains.nu = matrix(univt.mcmc$sims.array[,,'nu'], ncol=1)
chains.sigma = matrix(univt.mcmc$sims.array[,,'sigma'], ncol=1)

par(mfrow=c(2,2))
ts.plot(chains.mu, xlab="iteration", ylab="", main="mu")
ts.plot(chains.sigma, xlab="iteration", ylab="", main="sigma")
ts.plot(chains.nu, xlab="iteration", ylab="", main="df")

par(mfrow=c(2,2))
acf(chains.mu, main="mu")
acf(chains.sigma, main="sigma")
acf(chains.nu, main="df")

library(timeDate)
skewness(chains.nu)
kurtosis(chains.nu)

par(mfrow=c(2,2))
hist(chains.mu, main="mu")
hist(chains.sigma, main="sigma")
hist(chains.nu, main="df")

par(mfrow=c(2,2))
plot(density(chains.mu), main="mu")
plot(density(chains.sigma), main="sigma")
plot(density(chains.nu), main="df")

kurtosis.mcmc = 3*(chains.nu - 2) / (chains.nu - 4)
kurtosis.mcmc[chains.nu <= 4] = Inf
quantile(kurtosis.mcmc, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
sum(kurtosis.mcmc == Inf) / length(kurtosis.mcmc)

# Bootstrap estimate for kurtosis
library(MASS)
sampleSize = length(ibm) * 0.75
iterations = 1000
kurtosis.boot = rep(0, iterations)
for(i in 1:iterations){
  samples = sample(ibm, sampleSize, replace=T)
  estDF = as.numeric(fitdistr(samples, densfun='t')$estimate['df'])
  kurtosis.boot[i] = if(estDF <= 4) Inf else 3*(estDF - 2)/(estDF - 4)
}
quantile(kurtosis.boot, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
sum(kurtosis.boot == Inf) / length(kurtosis.boot)

########   AR 1, GDP data    ########
data(Tbrate, package="Ecdat")
# r = the 91-day treasury bill rate
# y = the log of real GDP
# pi = the inflation rate
diffTbrate = diff(Tbrate)
y = as.numeric(diffTbrate[,2])
y = y - mean(y)
N = length(y)
gdpDat = list("y", "N")

ar1_bugs_code="
model{
  # AR1 model y[i] = phi * y[i - 1] + e
  for(i in 2:N){
    y[i] ~ dnorm(mu[i], tau)
	mu[i] <- y[i-1]*phi
  }
  
  # Prior distribution
  phi ~ dnorm(0, .00001)
  tau ~ dgamma(0.1, .0001)
  sigma <- 1/sqrt(tau)
}
"
ar1_bugs_filename = "ar1.bug"
ar1_bugs_file = file(ar1_bugs_filename)
writeLines(ar1_bugs_code, ar1_bugs_file)
close(ar1_bugs_file)

inits.ar1 = function(){list(phi=rnorm(1, 0, .3), tau=runif(1, 1, 10))}
ar1.mcmc = bugs(gdpDat, inits.ar1, model.file=ar1_bugs_filename,
  parameters=c("phi", "sigma"), n.chains=3, n.iter=2600,
  n.burnin=100, n.thin=1, bugs.seed=5460, codaPkg=F,
  bugs.directory="E:/windows/programs/WinBUGS14")
file.remove(ar1_bugs_filename)
print(ar1.mcmc, digits=3)
plot(ar1.mcmc)
arima(y, order=c(1, 0, 0))


######  MA 1, simulated data  ######
set.seed(5640)
N = 800
y = arima.sim(n = N, list(ma=-.5), sd=.4)
y = as.numeric(y)
q = 5
ma.simData = list("y", "N", "q")

ma1_bugs_code="
model{
  for(i in 1:(N+q)){w[i] ~ dnorm(0, tau)}
  mu[1] <- w[1] + M
  for(i in 2:N){mu[i] <- w[i] + theta*w[i-1]}
  # The precision value can influence the estimation
  #   result greatly, try with different value instead
  #   of 10000, say 100 or 1000000, the estimation result
  #   will be very different.
  for(i in 1:N){y[i] ~ dnorm(mu[i], 10000)}
  
  # Prior distribution
  theta ~ dnorm(-.5, 0.00001)
  tau ~ dgamma(0.01, 0.01)
  sigma <- 1/sqrt(tau)
  M ~ dnorm(0, 0.001)
  for(i in 1:q){ypred[i] <- w[N+i] + theta*w[N+i-1]}
}
"
ma1_bugs_filename = "ma1.bug"
ma1_bugs_file = file(ma1_bugs_filename)
writeLines(ma1_bugs_code, ma1_bugs_file)
close(ma1_bugs_file)

inits.ma1 = function(){list(theta=rnorm(1, -.5, 1), tau=runif(1, 5, 8))}
ma1.mcmc = bugs(ma.simData, inits.ma1, model.file=ma1_bugs_filename,
  parameters=c("theta", "sigma", "ypred"), n.chains=3, n.iter=40000,
  n.burnin=10000, n.thin=5, codaPkg=F, bugs.seed=5460,
  bugs.directory="E:/windows/programs/WinBUGS14")
file.remove(ma1_bugs_filename)
print(ma1.mcmc, digits=3)
plot(ma1.mcmc)
par(mfrow=c(2,2))
acf(c(ma1.mcmc$sims.array[,,'theta']), main="theta")
acf(c(ma1.mcmc$sims.array[,,'sigma']), main="sigma")
acf(c(ma1.mcmc$sims.array[,,'ypred[1]']), main="ypred[1]")
acf(c(ma1.mcmc$sims.array[,,'ypred[2]']), main="ypred[2]")


######  ARMA(1,1), simulated data  ######
set.seed(5640)
N = 600
y = arima.sim(n=N, list(ar=.9, ma=-.5, sd=.4))
y = as.numeric(y)
arma.simData = list("y", "N")
ar1ma1_bugs_code="
model{
  for(i in 1:N){w[i] ~ dnorm(0, tau)}
  mu[1] <- w[1] + M
  for(i in 2:N){mu[i] <- phi * y[i - 1] + theta*w[i-1] + w[i]}
  # Note: the precision setup influences the estimation result
  for(i in 1:N){y[i] ~ dnorm(mu[i], 10000)}
  
  # Prior distribution
  phi ~ dnorm(0, .00001)
  theta ~ dnorm(-.5, 0.00001)
  tau ~ dgamma(0.01, 0.01)
  sigma <- 1/sqrt(tau)
  M ~ dnorm(0, 0.001)
}
"
ar1ma1_bugs_filename = "ar1ma1.bug"
ar1ma1_bugs_file = file(ar1ma1_bugs_filename)
writeLines(ar1ma1_bugs_code, ar1ma1_bugs_file)
close(ar1ma1_bugs_file)

inits.ar1ma1 = function(){list(phi=rnorm(1, 0, .3), theta=rnorm(1, -.5, 1), tau=runif(1, 5, 8))}
ar1ma1.mcmc = bugs(arma.simData, inits.ar1ma1, model.file=ar1ma1_bugs_filename,
  parameters=c("theta", "sigma", "phi"), n.chains=3, n.iter=40000,
  n.burnin=10000, n.thin=5, codaPkg=F, bugs.seed=5460,
  bugs.directory="E:/windows/programs/WinBUGS14")
file.remove(ar1ma1_bugs_filename)
print(ar1ma1.mcmc, digits=3)
arima(y, order=c(1,0,1))