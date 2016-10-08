library('Ecdat')
data(Mishkin)
# Plot inflation rate
plot(Mishkin[,2])
plot(diff(Mishkin[,2]))
plot(AirPassengers)
acf(Mishkin[,2], main="Inflation Rate")
acf(diff(Mishkin[,2]), main="Inflation Rate")
library('evir')
data(bmw)
acf(bmw)
Box.test(bmw, type="Ljung-Box",  lag=5)
Box.test(bmw, type="Box-Pierce", lag=5)
fitAR1 = arima(x=bmw, order=c(1, 0, 0))
resFitAR1 = residuals(fitAR1)
Box.test(resFitAR1, type="Box-Pierce", fitdf=-4)
acf(resFitAR1)
qqnorm(resFitAR1, datax=T)
plot(resFitAR1)
fitMishKin2 = arima(Mishkin[,2], order=c(1, 0, 0))
resFitMishKin2 = residuals(fitMishKin2)
Box.test(resFitMishKin2, type="L", lag=12)

library('forecast')
auto.arima(diff(Mishkin[,1]), max.p=20, max.q=0, max.P=0, max.Q=0, ic="aic")
auto.arima(diff(Mishkin[,1]), max.p=20, max.q=0, max.P=0, max.Q=0, ic="bic")
auto.arima(diff(Mishkin[,1]), max.p=0, max.q=20, max.P=0, max.Q=0, ic="aic")
auto.arima(diff(Mishkin[,1]), max.p=0, max.q=20, max.P=0, max.Q=0, ic="bic")

aicMishKin1 = rep(NA, 20)
bicMishKin1 = rep(NA, 20)
# auto regressive model, with first order differencial
arpFirstOrderModel = function(p){c(p, 1, 0)}
# moving average model, with first order differencial
maqFirstOrderModel = function(q){c(0, 1, q)}
arimaModelParam = maqFirstOrderModel
for(i in 1:20){
  aicMishKin1[i] =  arima(Mishkin[, 1], order=maqFirstOrderModel(i))$aic
  bicMishKin1[i] = aicMishKin1[i] + i * (log(length(diff(Mishkin[, 1]))) - 2)
}

plot(aicMishKin1, xlim=c(0,21), ylim=c(2450, 2550), pch='*', xlab="p", ylab="cretirian")
par(new=T)
plot(bicMishKin1, xlim=c(0,21), ylim=c(2450, 2550), pch='o', xlab="", ylab="")
legend("topright", c("AIC", "BIC"), pch=c('*', 'o'), inset=c(0.3, 0.0), bty="n")

library('Ecdat')
data(Capm)
estimates = matrix(data=NA, nrow=0, ncol=4, dimnames=list(NULL, c('p', 'q', 'aic', 'bic')))
for(p in 0:3){
  for(q in 0:3){
    aic =  arima(Capm[['rf']], order=c(p, 0, q))$aic
	aic = aic + 1290
    bic = aic + (p + q) * (log(length(Capm[['rf']])) - 2)
	estimates = rbind(estimates, c(p, q,  aic, bic))
  }
}
estimates
ar1MA1Est = arima(Capm[['rf']], order=c(1, 0, 1))
resAr1MA1Est = residuals(ar1MA1Est)
acf(resAr1MA1Est)
qqnorm(resAr1MA1Est, datax=T)
plot(resAr1MA1Est)


arimaSim = arima.sim(list(ar=c(0.5)), 500) # Simlate AR(1) model, with parameters: 0.5
arimaSim1stIntegral = cumsum(arimaSim)  # First order integral of the sample
arimaSim2ndIntegral = cumsum(arimaSim1) # Second order integral of the sample
par(mfrow=c(3,1))
plot(arimaSim, type="l")
plot(arimaSim1stIntegral, type="l")
plot(arimaSim2ndIntegral, type="l")

library('fEcofin')
data(CPI.dat)
cpiDates = as.Date(as.character(CPI.dat[[1]]))
indexRange = cpiDates >= as.Date('1977-01-01')
indexRange = indexRange & cpiDates <= as.Date('1987-12-31')
dateRange = as.Date(CPI.dat[[1]][indexRange])
logCPIRange = log(as.numeric(CPI.dat[[2]][indexRange]))
logCPIRange1stOrderDiff = diff(logCPIRange)
logCPIRange2ndOrderDiff = diff(logCPIRange1stOrderDiff)

par(mfrow=c(3,2))
plot(dateRange, logCPIRange, type="o", pch="o")
acf(logCPIRange)
plot(dateRange[2:length(dateRange)], logCPIRange1stOrderDiff, type="o", pch="o")
acf(logCPIRange1stOrderDiff)
plot(dateRange[3:length(dateRange)], logCPIRange2ndOrderDiff, type="o", pch="o")
acf(logCPIRange2ndOrderDiff)
logCPIAR0I2MA2 = arima(logCPIRange, order=c(0, 2, 2))
resLogCPIAR0I2MA2 = residuals(logCPIAR0I2MA2)
acf(resLogCPIAR0I2MA2)
Box.test(resLogCPIAR0I2MA2)

data(IP.dat)
par(mfrow=c(2, 2))
# Industrial Production
ip = IP.dat[[2]]
plot(ip)
ip1stOrderDiff = diff(ip)
plot(ip1stOrderDiff)
acf(ip1stOrderDiff)

library('forecast')
auto.arima(ip1stOrderDiff, max.p=20, max.q=0, max.P=0, max.Q=0, ic="aic")
ipEst = auto.arima(ip1stOrderDiff, max.p=3, max.q=3, d=0, max.P=0, max.Q=0, ic="aic")
acf(residuals(ipEst))

library(tseries)
adf.test(Mishkin[,1])
pred = predict(arima(Mishkin[,1], order=c(3, 1, 0)), n.ahead=100)$pred
allData = c(as.numeric(Mishkin[,1]), as.numeric(pred))

library(Ecdat)
library(tseries)
# r = the 91-day treasury bill rate
# y = the log of real GDP
# pi = the inflation rate
data(Tbrate)
plot(Tbrate)
acf(Tbrate)
adf.test(Tbrate[,1])
adf.test(Tbrate[,2])
adf.test(Tbrate[,3])

diff_rate = diff(Tbrate)
adf.test(diff_rate[,1])
adf.test(diff_rate[,2])
adf.test(diff_rate[,3])
pairs(diff_rate)
plot(diff_rate)
par(mfrow=c(1,1))
boxplot(diff_rate[,1] ~ cycle(diff_rate))

myData = Tbrate[,1]
# [AR, MA, Seasonal AR, Seasonal MA, Period, order of non-seasonal diff, order of seasonal diff]
fittedOrder = auto.arima(myData,max.P=0,max.Q=0,ic="aic")$arma
fit1 = arima(myData, order=c(fittedOrder[1], fittedOrder[6], fittedOrder[2]))
resFit1 = residuals(fit1)
acf(resFit1)
Box.test(resFit1, lag=10, type="L")
# GARCH effects, that is, volatility clustering, can be detected by looking for
# auto-correlation in the mean-centered squared residuals
resFit1Squared = resFit1^2
acf(resFit1Squared)
Box.test(resFit1Squared, lag=10, type="L")

myData = diff(Tbrate[,3])
fittedOrder = auto.arima(myData,max.P=0,max.Q=0,ic="aic")$arma
fit1 = arima(myData, order=c(fittedOrder[1], fittedOrder[6], fittedOrder[2]))
numPredict = 36
forecasts = predict(fit1, numPredict)
plot(myData, xlim=c(1980, 2006), ylim=c(-7, 12))
lines(seq(from=1997, by=.25, length=numPredict), forecasts$pred, col="red")
# Assuming Gaussian distribution on prediction error,
#  plotting the upper and lower limit, 0.95 percentile
lines(seq(from=1997, by=.25, length=numPredict), forecasts$pred + 1.96*forecasts$se, col="blue")
lines(seq(from=1997, by=.25, length=numPredict), forecasts$pred - 1.96*forecasts$se, col="blue")

library(Ecdat)
data(CRSPday)
crsp = CRSPday[,7]
par(mfrow=c(1,2))
acf(crsp)
acf(as.numeric(crsp))
arima(crsp, order=c(1,0,0))
arima(crsp, order=c(2,0,0))

data(Mishkin)
tb1 = log(Mishkin[,3])
par(mfrow=c(2,3))
acf(tb1, main=paste("diff=0"))
Box.test(tb1)
for(i in 1:5){
  tb1Diff = diff(tb1, differences=i)
  acf(tb1Diff, main=paste("diff=", as.character(i)))
  print(Box.test(tb1Diff))
}
aicFit = auto.arima(tb1, max.p=7, max.q=7, max.P=0, max.Q=0, ic="aic")
bicFit = auto.arima(tb1, max.p=7, max.q=7, max.P=0, max.Q=0, ic="bic")
aicFit
bicFit
Box.test(aicFit$residuals)
Box.test(bicFit$residuals)
par(mfrow=c(1,2))
acf(aicFit$residuals, main="acf for residuals using aic")
acf(bicFit$residuals, main="acf for residuals using bic")