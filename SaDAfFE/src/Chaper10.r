library(Ecdat)
plot(Hstarts)

hsts = Hstarts[,1]
acf(hsts)
boxplot(hsts~cycle(hsts))
par(mfrow=c(3,2))
plot(diff(hsts), type="o", pch="o", main="non-seasonal")
acf(diff(hsts), main="non-seasonal")
plot(diff(hsts, lag=frequency(hsts)), type="o", pch="o", main="seasonal")
acf(diff(hsts, lag=frequency(hsts)), main="seasonal")
plot(diff(diff(hsts, lag=frequency(hsts))), type="o", pch="o", main="seasonal&non-seasonal")
acf(diff(diff(hsts, lag=frequency(hsts))), main="seasonal&non-seasonal")
fit1 = arima(x=hsts, order=c(1,1,1), seasonal=list(order=c(1,1,1), period=4))
fit2 = arima(x=hsts, order=c(1,1,1), seasonal=list(order=c(0,1,1), period=4))
numPredict = 36
forecasts = predict(fit2, numPredict)

plot(hsts, xlim=c(1960, 2010), ylim=c(7.9, 12), type="o", pch="o")
lines(seq(from=2002, by=.25, length=numPredict), forecasts$pred, col="red", type="o", pch="o")
# Assuming Gaussian distribution on prediction error,
#  plotting the upper and lower limit, 0.95 percentile
lines(seq(from=2002, by=.25, length=numPredict), forecasts$pred + 1.96*forecasts$se, col="blue")
lines(seq(from=2002, by=.25, length=numPredict), forecasts$pred - 1.96*forecasts$se, col="blue")

library(fEcofin)
library(FitAR)
# Calculate the optimum transformation for the hoursing start series: Hstarts
BoxCox.Arima(fit1)
startDate = as.Date('1977-01-01')
endDate = as.Date('1988-01-01')
cpiDate = as.Date(CPI.dat[[1]])
ipDate = as.Date(IP.dat[[1]])
CPIDat1stOrderDiff = diff(log(CPI.dat[[2]][cpiDate > startDate & cpiDate < endDate]))
IPDat1stOrderDiff = diff(log(IP.dat[[2]][ipDate > startDate & ipDate < endDate]))
cpiVSip = as.data.frame(na.omit(cbind(CPIDat1stOrderDiff, IPDat1stOrderDiff)))
names(cpiVSip) = c("Delta CPI", "Delta IP")
ccf(CPIDat1stOrderDiff, IPDat1stOrderDiff)
arObj = ar(cpiVSip, order.max=1)
acf(na.omit(arObj$resid))

# TODO: write a multi-variate prediction function based on the arima estimate
arPredict = function(dat, mu, phi, num, sigma=NULL, numIter=1, confInterval = 0.95){
  columns = length(dat)
  if(columns != length(mu)) stop("number of mean values does not match number of variables")
  rows = length(dat[,1])
  numPhi = length(phi)
  predict_one = function(){
    pred = rep(mu, num)
    for(i in 1:num){
      k = (i - 1)*columns
      p = numPhi
      # pred = mu
	  startIndex = (i - 1) * columns + 1
	  endIndex = i * columns
      while(k > 0 && p > 0){
        pred[startIndex:endIndex] = pred[startIndex:endIndex] + phi[[p]] %*% (pred[(k - columns + 1):k] - mu)
        p = p - 1
        k = k - columns
      }
      k = rows
      while(k > 0 && p > 0){
        pred[startIndex:endIndex] = pred[startIndex:endIndex] + phi[[p]] %*% (as.numeric(dat[k,]) - mu)
        p = p - 1
        k = k - 1
      }
	  if (!is.null(sigma)) pred[startIndex:endIndex] = mvrnorm(n=1,mu=pred[startIndex:endIndex],Sigma=sigma)
    }
    pred
  }
  if(is.null(sigma)) numIter = 1
  predicted = matrix(nrow = numIter, ncol=columns * num)
  for(j in 1:numIter){
      predicted[j,] = predict_one()
  }
  lower = NULL
  upper = NULL
  if(confInterval >= 1e-6 && confInterval <= (1 - 1e-6)){
    lowQantile = function(x) {
	  quantile(x, (1 - confInterval) / 2)
	}
	upQantile = function(x) {
	  quantile(x, (1 + confInterval) / 2)
	}
	lower = matrix(apply(predicted, 2, lowQantile), nrow=num, byrow=T)
	upper = matrix(apply(predicted, 2, upQantile), nrow=num, byrow=T)
  }
  predicted = matrix(apply(predicted, 2, mean), nrow=num, byrow=T)
  list(pred=predicted, l=lower, u=upper)
}
phi = arObj$ar[,,]
predicted = arPredict(cpiVSip, arObj$x.mean, list(phi), 20, arObj$var.pred, 50000)

plot(1:20,predicted$pred[,1],type="b",lty=2,ylim=c(0,.006),pch="*",cex=3,
   xlab="k",ylab="forecast",cex.axis=1.5,cex.lab=1.5,lwd=3)

lines(1:20,predicted$pred[,2],type="b",pch="o",cex=2,lwd=3)
abline(h=arObj$x.mean[1])
abline(h=arObj$x.mean[2],lty=2)
legend(2,.0015,c("CPI","IP"),lty=c(2,1),pch=c("*","o"),cex=1.6,lwd=3,
  pt.cex=c(3,2))

par(mfrow=c(2,1))
plot(1:20,predicted$u[,1],type="b",lty=2,ylim=c(-0.003,.013),lwd=2,
   xlab="lag",ylab="forecast",cex.axis=1.5,cex.lab=1.5,
   main="CPI",cex=1.5)
lines(1:20,predicted$l[,1],type="b",lty=2,lwd=2,cex=1.5)
lines(1:20,predicted$pred[,1],type="b",lwd=2,pch="*",cex=1.5)
plot(1:20,predicted$u[,2],,type="b",lty=2,ylim=c(-0.02,.025),lwd=2,
   xlab="lag",ylab="forecast",cex.axis=1.5,cex.lab=1.5,
   main="IP",cex=1.5)
lines(1:20,predicted$l[,2],type="b",lty=2,lwd=2,pch="o",cex=1.5)
lines(1:20,predicted$pred[,2],type="b",lwd=2,pch="*",cex=1.5)


library(longmemo)

D = c(-.35, .35, 0.7)
set.seed("09201948")
par(mfrow=c(length(D),2))
for (i in 1:length(D))
{
  H = D[i] + 0.5
  A = as.integer(H)
  H = H - A
  x = simARMA0(2500,H)
  while(A > 0){
    x = cumsum(x)
  A = A - 1
  }
  plot(x,main=toString(paste("d =", D[i])), type="l")
  acf(x,main=toString(paste("d =", D[i])))
}

y = simARMA0(2500,.7 - 1 + 1/2)
x = cumsum(y)
par(mfrow=c(2,2))
plot(x,main="d = 0.7",type="l")
acf(x,main="d = 0.7")

library(fracdiff)
acf(diffseries(x,.7),main = expression(paste(Delta^0.7, Y)))
acf(diff(x),main = expression(paste(Delta, Y)))

par(mfrow=c(2,2))
acf(Mishkin[,1])
acf(diff(Mishkin[,1]))
acf(diffseries(Mishkin[,1], 0.4))

library(Ecdat)
data(IncomeUK)
consumption = IncomeUK[,2]
plot(consumption)
fitConsumption = auto.arima(consumption)
# > fitConsumption
# Series: consumption 
# ARIMA(0,1,0)(1,1,0)[4]                    
# 
# Coefficients:
#          sar1
#       -0.3904
# s.e.   0.1308
# 
# sigma^2 estimated as 212648:  log likelihood=-400.11
# AIC=804.23   AICc=804.47   BIC=808.17
acf(fitConsumption$residuals)
logConsumption = log(consumption)
fitLogConsumption = auto.arima(logConsumption)
# > fitLogConsumption
# Series: log(consumption) 
# ARIMA(0,1,0)(0,1,1)[4]                    
# 
# Coefficients:
#          sma1
#       -0.5348
# s.e.   0.1164
# 
# sigma^2 estimated as 0.0003056:  log likelihood=139.77
# AIC=-275.55   AICc=-275.31   BIC=-271.61
acf(fitLogConsumption$residuals)
fitLogConsumptionBIC = auto.arima(logConsumption, ic="bic")
# > fitLogConsumptionBIC
# Series: log(consumption) 
# ARIMA(0,1,0)(0,1,1)[4]                    
# 
# Coefficients:
#          sma1
#       -0.5348
# s.e.   0.1164
# 
# sigma^2 estimated as 0.0003056:  log likelihood=139.77
# AIC=-275.55   AICc=-275.31   BIC=-271.61
acf(fitLogConsumptionBIC$residuals)
forecastLogConsumption = forecast(fitLogConsumptionBIC, h=8)

data(Tbrate, package="Ecdat")
# r = the 91-day Treasury bill rate
# y = the log of real GDP
# pi = the inflation rate
rate1stOrderDiff = diff(Tbrate)
arFitTRate = ar(rate1stOrderDiff, order.max=4, aic=T)
# > arFitRate
# 
# Call:
# ar(x = rate1stOrderDiff, aic = T, order.max = 4)
# 
# $ar
# , , 1
# 
#             r       y         pi
# r   0.2139366 17.2841 -0.0304881
# y  -0.0001014  0.2206  0.0001002
# pi  0.1846333 16.1098 -0.2487303
# 
# 
# $var.pred
#           r         y        pi
# r  0.807888  0.001858  0.150369
# y  0.001858  0.000138 -0.001082
# pi 0.150369 -0.001082  3.682619
acf(arFitTRate$resid)
phi = arFitTRate$ar[,,]
predicted = arPredict(as.data.frame(rate1stOrderDiff), arFitTRate$x.mean, list(phi), 20)

data(Mishkin, package="Ecdat")
cpi = as.vector(Mishkin[, 5])
sqrtCPI1stOrder = diff(sqrt(cpi))
acf(sqrtCPI1stOrder)

library('fracdiff')
fit.frac = fracdiff(sqrtCPI1stOrder, nar=0, nma=0)
fit.frac$d
fdiff = diffseries(sqrtCPI1stOrder, fit.frac$d)
acf(fdiff)

library(forecast)
arFDiff = auto.arima(fdiff)
acf(arFDiff$residuals)

library(AER)
library(forecast)
data("FrozenJuice")
price = FrozenJuice[,1]
plot(price)
fit = auto.arima(price, ic="bic")

phi = fit$model$phi
n = length(price)
m = length(phi)
sink("priceBootstrap.txt")
set.seed(1998852)
# > fit
# Series: price 
# ARIMA(2,1,0)                    
# 
# Coefficients:
#          ar1     ar2
#      0.2825  0.0570
# s.e.  0.0407  0.0408
# 
# sigma^2 estimated as 10.02:  log likelihood=-1570.11
# AIC=3146.23   AICc=3146.27   BIC=3159.47
for (iter in 1:10) {
  eps = rnorm(n + 20)
  y = rep(0, n + 20)
  for (t in (m + 1) : (n + 20)){
    # phi = 0.2825  0.0570
    y[t] =  phi %*% y[(t - m) : (t - 1)] + eps[t]
  }
  # order is 1
  y = ts(cumsum(y[101:n + 20]), frequency=12)
  fit = auto.arima(y, d=1, D=0, ic="bic")
  print(fit)
}
sink()

set.seed(1998852)
niter = 250
estimates = matrix(0, nrow=niter, ncol=m)
for(iter in 1:niter) {
  eps = rnorm(n + 20)
  y = rep(0, n + 20)
  for (t in (m + 1):(n + 20)) {
    y[t] = phi %*% y[(t - m) : (t - 1)] + eps[t]
  }
  y = ts(cumsum(y[101:n + 20]), frequency=12)
  fit = arima(y, order=c(2, 1, 0))
  estimates[iter,]=fit$coef
}

apply(estimates, 2, mean)
apply(estimates, 2, var)

library(Ecdat)
data(IncomeUK)
income = IncomeUK[,1]
nonSeasonal = auto.arima(income, max.P=0, max.Q=0, max.D=0)
seasonal = auto.arima(income)
nonSeasonal
seasonal
par(mfrow=c(1, 2))
acf(nonSeasonal$residuals)
acf(seasonal$residuals)

data(USMacroG, package="AER")
unemployment = USMacroG[, "unemp"]
fitUnemp = auto.arima(unemployment)
plot(unemployment, ylim=c(3, 12))
par(new=T)
plot(unemployment - fitUnemp$residuals, col="red", ylim=c(3, 12))
# > fitUnemp
# Series: unemployment 
# ARIMA(1,1,0)(0,0,2)[4]                    
# 
# Coefficients:
#         ar1     sma1     sma2
#      0.6611  -0.4199  -0.2623
# s.e.  0.0544   0.0675   0.0660
# 
# sigma^2 estimated as 0.08112:  log likelihood=-32.85
# AIC=73.69   AICc=73.89   BIC=86.94
# > fitUnemp$sigma2
# [1] 0.08112103

# instead of writing general code, use the output of the above to test
# the correctness of the estimator

n = length(unemployment)
for (iter in 1:10) {
  # standard deviation of error
  sdE = sqrt(fitUnemp$sigma2)
  eps = rnorm(n + 20, sd=sdE)
  y = rep(0, n + 20)
  # non-seasonal
  # non-seasonal AR parameters
  phi = c(0.6611)
  # p is number of parameters for non-seasonal ARIMA
  ############################
  #  data from the non-seasonal ARIMA(1,1,0) process
  p = length(phi)
  for (t in (p + 1) : (n + 20)){
    y[t] =  phi %*% y[(t - p) : (t - 1)] + eps[t]
  }
  # order is 1
  # frequency is 4
  y = ts(cumsum(y), frequency=4)
  ############################
  #  data from the seansonal ARIMA(0,0,2) process
  #  Since p for autoregressive is zero, we only have random noise here.
  sma1 = rnorm(n + 20, sd=sdE)
  sma2 = rnorm(n + 20, sd=sdE)
  y = y + (-0.4199) * sma1 + (-0.2623) * sma2
  fit = auto.arima(y[10:n+20], ic="bic")
  print(fit)
  print(fit$sigma2)
}