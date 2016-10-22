X11()

data(Mishkin, package="Ecdat")
data(SP500, package="Ecdat")
data(Garch, package="Ecdat")
data(Capm, package="Ecdat")
genYears = function(start, end, datas){
  n = if(is.null(dim(datas))) length(datas) else dim(datas)[1]; start + (1:n) * (end - start)/n
}

yearsData = list(
  Capm=genYears(1960, 2003, Capm), SP500=genYears(1981, 1991.25, SP500),
  Mishkin=genYears(1950+1/12, 1991, Mishkin), Garch = genYears(1980, 1987.5, Garch))

garchPlot = function(colName, target){
  target[[paste("diff", colName, sep="")]] = list(dat=as.vector(diff(Garch[, colName])), 
    years = (yearsData$Garch)[-1], params=list(ylab="|change in rate|",
	  main=paste(colName, "/dollar exchange rate", sep="")))
  target
}
genPlot = function(dat, yearsCol, params){
  list(dat=dat, years=yearsData[[yearsCol]], params=params)
}

garchPlots = function(colNames){
  plots = as.list(NULL)
  for(i in 1:length(colNames)){plots = garchPlot(colNames[i],plots)}
  plots
}

plots = garchPlots(c("dm", "bp", "cd", "dy"))
plots$difflogrf = genPlot(diff(log(Capm$rf)), "Capm", 
  list(main="Risk-free interest rate", ylab="|change in log(rate)|"))
plots$difflogrf$years = plots$difflogrf$years[-1]
plots$infl = genPlot(scale(as.vector(Mishkin[,1])), "Mishkin",
  list(ylab="|rate - mean(rate)|", main="Inflation rate", span=.3))
plots$r500= genPlot(SP500$r500, "SP500",
  list(ylab="|log return|", main="S&P 500 daily return"))

par(mfrow=c(2,4))
for(i in 1:length(plots)){
  plt = plots[[i]]; dat = abs(plt$dat); param = plt$params
  plot(plt$years, dat, main=param$main, type="l", xlab="years",
    cex.axis=1.5, cex.lab=1.5, cex.main=1.5, ylab=param$ylab, ylim=param$ylim)
  mod = loess(dat~plt$years,span=if(is.null(param$span)) .25 else param$span)
  lines(plt$years, predict(mod), lwd=6, col=gray(.7))
}

plot_all = function(e, a, u, g=F){
  index = 0
  getMain = function(gt=if(g) "GARCH(1,1)" else "ARCH(1)", p="", ar=F, sqrt=F){
    paste("(", letters[index <<- index + 1], ") ", p, if(ar) "AR(1)/" else "",
	  gt, if(sqrt) " squared" else "", sep="")
  }
  plot(e[11:n], type="l", xlab="t", ylab=expression(epsilon), main=getMain("white noise"))
  plot(sqrt(sig2[11:n]), type="l", xlab="t", ylab=expression(sigma[t]), 
    main=getMain(p="conditional std dev"))
  plot(a[11:n], type="l", xlab="t", ylab="a", main=getMain())
  plot(u[11:n], type="l", xlab="t", ylab="u", main=getMain(ar=T))
  acf(a[11:n], main=getMain())
  acf(a[11:n]^2, main=getMain(sqrt=T))
  acf(u[11:n], main=getMain(ar=T))
  acf(u[11:n]^2, main=getMain(ar=T, sqrt=T))
}

#iterations = 1000
#variations = rep(0, iterations)
#for(j in 1:iterations){
n=10010; e=rnorm(n); a=e; u=e; sig2=e^2; omega=1; alpha1=.95; ar1=.8; mu=.1

par(mfrow=c(4,4))
# a is ARCH(1) process
# u is AR(1)/ARCH(1) process
for(i in 2:n){
  sig2[i] = omega + alpha1 * a[i - 1]^2
  a[i] = sqrt(sig2[i]) * e[i]
  u[i] = mu + ar1*(u[i - 1] - mu) + a[i]
}
plot_all(e, a, u)

alpha1 = .08; beta1 = .9 #; ma1 = 
# a is GARCH(1, 1) process
# u is AR(1)/GARCH(1, 1) process
for(i in 2:n){
  sig2[i] = omega + alpha1 * a[i - 1]^2 + beta1*sig2[i - 1]
  a[i] = sqrt(sig2[i]) * e[i]
  u[i] = mu + ar1*(u[i - 1] - mu) + a[i] # + ma1 * a[i - 1]
}
plot_all(e, a, u, T)

#variations[j] = var(a)
#}

# by theory, the unconditional variance of a should be : omega/(1 - alpha1) when alpha1 < 1 ?
# with 1000 iterations of above loop, it seems var(a) != omega/(1 - alpha1) when alpha1 < 1
# maybe it's because the variability is too large.
var(a)

library(fGarch)
library(MASS)
data(bmw, package="evir")

garchNorm = garchFit(~arma(1, 0) + garch(1, 1), data=bmw, cond.dist="norm")
garchT = garchFit(~arma(1, 1) + garch(1, 1), data=bmw, cond.dist="std")
apgarchT = garchFit(~arma(1,0) + aparch(1,1), include.delta=T, cond.dist="std",data=bmw)
qqResid = function(model){
  summary(model)
  stdResid = model@residuals / model@sigma.t
  callStr = paste(model@call$formula[2], model@call$cond.dist, "residual")
  qqnorm(stdResid, datax=T, ylab="Standardized residual quantiles",
    main=paste("(a)", callStr, "normal plot"), xlab="normal quantiles")
  qqline(stdResid, datax=T)
  # tVars = rt(length(stdResid), df=4)
  tVars = qt((1:length(stdResid))/(length(stdResid) + 1), df=4)
  qqplot(stdResid, tVars, main=paste("(b)", callStr, "t plot, df=4"), 
    xlab="Standardized residual quantiles", ylab="t quantiles")
  abline(lm(qt(c(.25, .75), df=4)~quantile(stdResid, c(.25, .75))))
  acf(stdResid, main=callStr)
}

fitdistr(garchNorm@residuals / garchNorm@sigma.t, densfun="t")

options(digits=10)
par(mfrow=c(3,3))
qqResid(garchNorm)
qqResid(garchT)
qqResid(apgarchT)

genData = function(dat, target=.5){
  rho = function(a, b, interval=0:9){
    a * (1 - a * b - b^2) / (1 - 2 * a * b - b^2) * (a + b)^interval
  }  
  result = list(NULL)
  for(i in 1:length(dat)){
    beta = uniroot(f=function(b){rho(dat[i], b, 0) - target}, interval=c(0, 1 - dat[i]))$root
    result[[paste("label", i, sep="")]] =
	  parse(text=paste('paste(alpha,"=",', dat[i], ',"," ,beta,"=",', beta, ')'))
	result[[paste("dat", i, sep="")]] = c(1, rho(dat[i], beta))
  }
  result
}
plotData = genData(c(.1, .3, .5))

plot(0:10, plotData$dat1, type="b", ylim=c(0, 1), lty=1, lwd=2, ylab=expression(paste(rho[a^2], "(lag)")))
lines(0:10, plotData$dat2, type="b", lty=2, lwd=2)
lines(0:10, plotData$dat3, type="b", lty=3, lwd=2)
legend("topright", c(plotData$label1, plotData$label2, plotData$label3), lty=1:3, lwd=2)


data(bmw,package="evir")
acf(residuals(arima(bmw,order=c(1,0,0)))^2)

gammas = c(-0.5, -0.2, 0, 0.12, 0.3, 0.9)
par(mfrow=c(3, 2))
for(i in 1:length(gammas)){
  x=c(-3, 0, 3)
  plot(x, abs(x) - gammas[i]*x, type='l', lwd=2, 
    main=parse(text=paste('paste(gamma,"=",', gammas[i], ')')),
    ylab=expression(paste(g[gamma], "(x)")))
}

library("fGarch")
library("Ecdat")
library("fEcofin")
library("forecast")
data(nelsonplosser)
names(nelsonplosser)
newData = na.omit(nelsonplosser)
newData = list(difflogsp=diff(log(newData$sp)), diffgnpr=diff(newData$gnp.r),
  diffcpi=diff(newData$cpi), diffbnd=diff(newData$bnd))
newData = as.data.frame(newData)

fitLm2 = lm(difflogsp ~ diffgnpr + diffcpi + diffbnd, data=newData)
resLm2 = rstudent(fitLm2)
auto.arima(resLm2)
xregressors = cbind(newData$diffgnpr, newData$diffcpi, newData$diffbnd)/100
fitArma = arima(newData$difflogsp, xreg=xregressors, order=c(0, 0, 1))
par(mfrow=c(1, 2))
acf(fitArma$resid)
acf(fitArma$resid^2)

nelploss.garch = garchFit(~arma(0, 1)+garch(1, 0), data=residuals(fitLm2))
summary(nelploss.garch)
nelploss.garch.std.resid = nelploss.garch@residuals / nelploss.garch@sigma.t
qqnorm(nelploss.garch.std.resid, datax=T)
qqline(nelploss.garch.std.resid, datax=T)

par(mfrow=c(2,2))
acf(resLm2, main="(a) regression: residuals")
acf(resLm2^2, main="(b) regression: squared residuals")
acf(nelploss.garch.std.resid, main="(c) MA/ARCH: residuals")
acf(nelploss.garch.std.resid^2, main="(d) MA/ARCH: squared residuals")

fitLm3 = lm(difflogsp ~ diffgnpr + diffcpi + diffbnd, data=newData,
  weights=1/(nelploss.garch@sigma.t^2))
summary(fitLm3)
plot(fitted(fitLm2), fitted(fitLm3))

data(bmw, package="evir")
year = genYears(1973, 1996+7/12, bmw)
origin1 = 4100; origin2 = 3800
myGarchFit = function(data){garchFit(~arma(1, 0) + garch(1, 1), data=data, cond.dist="std")}
myPredict = function(fit, alpha=0.05, nahead=1500){
  pred = predict(fit, n.ahead=nahead)
  t1 = qstd(1-alpha/2, mean=0, sd=1, nu=coef(fit)["shape"])
  list(pred=pred, ul=pred$meanForecast + t1*pred$standardDeviation,
    ll=pred$meanForecast - t1*pred$standardDeviation)
}
bmw.garchT1 = myGarchFit(bmw[1:origin1])
bmw.garchT1.pred = myPredict(bmw.garchT1)
bmw.garchT2 = myGarchFit(bmw[1:origin2])
bmw.garchT2.pred = myPredict(bmw.garchT2)

plotPred = function(years, pred, lwd=4, ...){
  lines(years, pred$ul, col="black", lwd=lwd, ...)
  lines(years, pred$ll, col="black", lwd=lwd, ...)
}

plot(year, bmw, type="l", xlab="year", ylab="return",
  xlim=c(1986, 1992), ylim=c(-.13, .21), main="Forecasting BMW returns")
plotPred(year[(origin1+1):(origin1+nahead)], bmw.garchT1.pred)
plotPred(year[(origin2+1):(origin2+nahead)], bmw.garchT2.pred, lty=2)
legend("topleft", c("11-15-87", "9-18-88"), lty=c(2, 1), lwd=4, bty="n")

data(Tbrate, package="Ecdat")
library(tseries)
library(fGarch)

# r = the 91-day treasury bill rate
# y = the log of real GDP
# pi = the inflation rate
Tbill = Tbrate[, 1]
deltaTbill = diff(Tbill)
deltaLogTbill = diff(log(Tbill))

par(mfrow=c(2,2))
plot(Tbill)
acf(Tbill)
plot(deltaTbill)
acf(deltaTbill)
adf.test(Tbill)
adf.test(deltaTbill)
kpss.test(Tbill)
kpss.test(deltaTbill)

garchDeltaTbill = garchFit(formula= ~arma(1, 0) + garch(1, 1), deltaTbill)
garchTbill = garchFit(formula=~arma(1, 0) + garch(1, 1), Tbill)
garchDeltaLogTbill = garchFit(formula=~arma(1, 0) + garch(1, 1), deltaLogTbill)

myPlot = function(fit, ...){
  res = residuals(fit)
  resStd = res / fit@sigma.t
  plot(res, ...)
  acf(res, ...)
  acf(res^2, ...)
  plot(resStd, ...)
  acf(resStd, ...)
  acf(resStd^2, ...)
}

par(mfrow=c(3, 6))
myPlot(garchDeltaTbill, main="AR(1)/GARCH(1, 1) delta Tbill")
myPlot(garchTbill, main="AR(1)/GARCH(1, 1) Tbill")
myPlot(garchDeltaLogTbill, main="AR(1)/GARCH(1, 1) delta log Tbill")

library(fGarch)
data(SP500, package="Ecdat")

# index for the minum return(black Monday)
returnBlackMonday = min(SP500$r500)
blackMonIndex = (1:length(SP500$r500))[SP500$r500 == returnBlackMonday]
tradingDaysPerYear = 253
return2YearsBeforeBlackMonday = SP500$r500[(blackMonIndex - 2 * tradingDaysPerYear):(blackMonIndex - 1)]
plot(c(return2YearsBeforeBlackMonday, returnBlackMonday))
#results = garchFit(~arma(1, 0) + garch(1, 1), data=return2YearsBeforeBlackMonday,
results = garchFit(~arma(1, 0) + aparch(1, 1), data=return2YearsBeforeBlackMonday,
  cond.dist="std")
dfhat = as.numeric(results@fit$par['shape'])
forecast = predict(results,n.ahead=1)
summary(results)

pstd(returnBlackMonday, mean=as.numeric(forecast['meanForecast']),
  sd=as.numeric(forecast['standardDeviation']), nu=dfhat)
stdRes = residuals(results) / results@sigma.t
arFit = arima(return2YearsBeforeBlackMonday, order=c(1, 0, 0))
res = residuals(arFit)
par(mfrow=c(2, 3))
plot(stdRes, type='l')
acf(stdRes)
qqnorm(stdRes)
qqline(stdRes)

plot(res)
acf(res)
qqnorm(res)
qqline(res)

data(Irates, package='Ecdat')
diffr = as.numeric(diff(log(Irates[, 2])))
garchDiffRate = garchFit(~arma(1,0)+garch(1,1),data=diffr, cond.dist="std")
summary(garchDiffRate)

# Because: alpha1 + beta1 > 1, the variance of residual squared is not convergent.
#  so there's no formula for the acf of the residual squared
alpha1 = garchDiffRate@fit$par['alpha1']
beta1 = garchDiffRate@fit$par['beta1']

# With this fit alpha1 + beta1 < 1.
garchDiffRate = garchFit(~arma(1, 1) + garch(1, 1),data=diffr, cond.dist="sged")
param = garchDiffRate@fit$par
param
# >            mu           ar1           ma1         omega        alpha1         beta1          skew         shape 
# >  0.0083149348 -0.5181016314  0.5863487451  0.0003905782  0.2553392134  0.7249090484  0.8825486869  1.0000000000

stdResid = residuals(garchDiffRate)/garchDiffRate@sigma.t
sgedVars = rsged((1:length(stdResid))/(length(stdResid) + 1), nu=param['shape'], xi=param['skew'])
qqplot(stdResid, sgedVars, xlab="Standardized residual quantiles", ylab="sged quantiles")
abline(lm(qsged(c(.25, .75), nu=param['shape'], xi=param['skew'])~quantile(stdResid, c(.25, .75))))

library('mnormt')