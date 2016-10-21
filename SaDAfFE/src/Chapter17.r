datNoOmit = read.table("F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 17\\treasury_yields.txt",header=T)
diffdatNoOmit = diff(as.matrix(datNoOmit[,2:12]))
dat = na.omit(datNoOmit)
diffdat = na.omit(diffdatNoOmit)

n = dim(diffdat)[1]
options(digits=5)
pca = prcomp(diffdat)

summary(pca)

# eigenValVec = eigen(var(diffdat))
# rotations = diffdat %*% eigenValVec$vectors

time = c(1/12,.25,.5,1, 2, 3, 5, 7, 10, 20, 30)
par(mfrow=c(2,2))
my_plot1 = function(i1, i2, i3, xlab, ylab, main, ylim, ...){
  labls = c(my_plot(plot, i1, 1, ylim=ylim, ylab=ylab,xlab=xlab,main=main, ...), my_plot(lines, i2, 2), my_plot(lines, i3, 3))
  legend("bottomright", labls, lty=c(1,2,3),lwd=2)
}

par(mfrow=c(2,2))
my_plot = function(func, index, lty, ...){
  func(time, as.vector(dat[index, 2:(length(time)+1)]), type="b", lty=lty, lwd=2, ...)
  as.character(dat[index, 1])
}
my_plot1(1, 486, n+2, "T", "Yield", "(a)", c(0,6))
plot(pca,main="(b)")

my_plot = function(func, index, lty, ...){
  func(time, pca$rotation[,index], type="b", lty=lty, lwd=2, ...)
  paste("PC ", index)
}
my_plot2 = function(main, ...){
  my_plot1(1, 2, 3, "T", "Yield", main, c(-.8, .8), ...)
  lines(0:30,0*(0:30),lwd=1)
}

my_plot2("(c)")
my_plot2("(d)", xlim=c(0,3))

mu = apply(dat[, 2:(length(time)+1)], 2, mean)
plot_mean_pc = function(index, main, ...){
  plot(time, mu, ylim=c(0,6), type="b", lwd=4, xlab="T", ylab="Yield", main=main, ...)
  lines(time, mu+pca$rotation[,index], lty=5, type="b", lwd=2)
  lines(time, mu-pca$rotation[,index], lty=3, type="b", lwd=2)
  labls = c("mean", paste("mean + PC", index), paste("mean -  PC", index))
  legend("bottomright", labls, lty=c(1,5,3), lwd=c(4,2,2))
}

par(mfrow=c(2,2))
plot_mean_pc(1, "(a)", xlim=c(0, 7))
plot_mean_pc(2, "(b)", xlim=c(0, 7))
plot_mean_pc(3, "(c)", xlim=c(0, 7))

plot(time, pca$rotation[,4], ylim=c(-1,1), type="b", lwd=2, ylab="PC", xlab="T", xlim=c(0,30), main="(d)")
lines(time,pca$rotation[,5], lty=2, type="b", lwd=2)
lines(0:30, 0*(0:30), lwd=1)
legend("topright", c("PC 4","PC 5"), lty=c(1,2), lwd=2)

par(mfrow=c(1,3))
for (i in 1:3){
  plot(pca$x[,i], main=paste("PC", i), xlab="day", ylab="")
}

acf(pca$x[,1:3], ylab="", xlab="lag")

library(fEcofin)
startDate = as.Date('1977-06-01')
endDate = as.Date('1988-01-01')
cpiDate = as.Date(CPI.dat[[1]])
ipDate = as.Date(IP.dat[[1]])
CPI = as.data.frame(diff(log(CPI.dat[[2]][cpiDate > startDate & cpiDate < endDate])))
IP = as.data.frame(diff(log(IP.dat[[2]][ipDate > startDate & ipDate < endDate])))
names(CPI) = "CPI"
names(IP) = "IP"
arFit = ar(cbind(CPI, IP))

res = arFit$resid[6:125,]
berndt = as.matrix(berndtInvest[,3:11])
lmfit = lm(berndt~res[,1]+res[,2])
slmfit = summary(lmfit)
rsq = rep(0, ncol(berndt))
for(i in 1:ncol(berndt)){rsq[i] = slmfit[[i]][[8]]}
beta_CPI = lmfit$coef[2, ]
beta_IP = lmfit$coef[3, ]

par(mfrow=c(1, 3))
barplot(rsq, hori=T, names=names(beta_CPI), main="R squared")
barplot(beta_CPI, hori=T, main="beta CPI")
barplot(beta_IP, hori=T, main="beta IP")

data(CRSPmon, package="Ecdat")
FF_data = read.table("F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 17\\FamaFrench_mon_69_98.txt",header=T)
stockNames = c("ge", "ibm", "mobil")
famaFrenchFactorNames = c("Mkt.RF", "SMB", "HML")
stocks = as.matrix(CRSPmon[, stockNames]*100 - FF_data[, "RF"])
famaFrenchFactors = FF_data[, famaFrenchFactorNames]
allData = cbind(stocks, as.matrix(famaFrenchFactors))
colnames(stocks) = stockNames
colnames(allData) = c(stockNames, famaFrenchFactorNames)
pairs(allData)

fit = lm(stocks~Mkt.RF+SMB+HML, data=famaFrenchFactors)
fit
summary(fit)

cor(fit$resid)
cor.test(fit$resid[,1], fit$resid[, 2])
cor.test(fit$resid[,1], fit$resid[, 3])
cor.test(fit$resid[,2], fit$resid[, 3])

pairs(fit$resid)
library(robust)
covRob(fit$resid, corr=T)
# covRob(cbind(fit$resid, as.matrix(famaFrenchFactors)), corr=T)

sigmaFactor = var(famaFrenchFactors)
beta0 = fit$coef[1, ]
betaFactor = fit$coef[-1, ]
varError = var(fit$resid)
sigmaError = diag(diag(varError))
estimateCov1 = t(betaFactor) %*% sigmaFactor %*% betaFactor
estimateCov = estimateCov1 + sigmaError

library("fEcofin")
codes = as.matrix(model.matrix(~as.factor(c(3,3,2,1,1,2, 3,3,1,2,2,3,1,2,3))))
codes[, 1] =  1 - codes[, 2] - codes[, 3]
codes[, 2:3] = codes[, 1:2]
codes[, 1] = 1
returns = berndtInvest[, -c(1, 11, 18)]
betas = as.data.frame(codes[1:15, 1:3])
colnames(betas) = c("intercept", "tech", "oil")

betas

factors = matrix(0,nrow=120,ncol=3)
for (i in 1:120){
  return_data = cbind(t(returns[i, ]), betas)
  colnames(return_data)[1] = "return"
  lmfit = lm(return~.-1, data=return_data)
  factors[i, ]=lmfit$coef
}

print(sd(factors), digits=4)
par(mfrow=c(1, 3), cex.axis=1.08, cex.lab=1.08, cex.main=1.05)
my_plot = function(i, m){plot(factors[, i], main=m, lty="dotted", type="b", lwd=2, xlab="month", ylab="factor")}
my_plot(1, "market")
my_plot(2, "technology")
my_plot(3, "oil")

colnames(factors) = c("market", "tech", "oil")
cor(factors)
sqrt(diag(cov(factors)))

acf(factors, ylab="", xlab="lag")

data(equityFunds, package="fEcofin")
options(digits=3)
fa = factanal(equityFunds[, 2:9], 4, rotation="none")
fa

betas = fa$loadings[,]
estCor = betas %*% t(betas) + diag(fa$unique)
corDiff = estCor - fa$corr
max(corDiff)
min(corDiff)
eigDiff = eigen(corDiff)
sort(eigDiff$values)
wghts = matrix(1/8,nrow =1,ncol=8)
wghts %*% estCor %*% t(wghts)
wghts %*% fa$corr %*% t(wghts)


fa1 = factanal(equityFunds[,2:9],4,rotation="varimax")
fa1

wghts %*% (fa1$loadings[,] %*% t(fa1$loadings[,]) + diag(fa1$unique)) %*% t(wghts)

yieldDat = read.table("F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 17\\yields.txt", header=T)
maturity = c((0:5), 5.5, 6.5, 7.5, 8.5, 9.5)
pairs(yieldDat)
par(mfrow=c(4, 3))
for(i in 0:11){plot(maturity, yieldDat[100*i + 1, ], type="b")}
eig = eigen(cov(yieldDat))
eig$values
eig$vectors
par(mfrow=c(1,1))
barplot(eig$values)

par(mfrow=c(2,2))
for(i in 1:4){
  plot(eig$vectors[, i], ylim=c(-.7, .7), type="b", ylab=paste("eigen vector", i))
  abline(h=0)
}

par(mfrow=c(4,3))
for(i in 1:ncol(yieldDat)){
  plot(yieldDat[, i], type='l', ylab = names(yieldDat)[i])
}

library(tseries)
adf.test(yieldDat[, 1])

n = dim(yieldDat)[1]
diffYield = yieldDat[-1, ] - yieldDat[-n, ]
par(mfrow=c(4,3))
for(i in 1:ncol(yieldDat)){
  plot(diffYield[, i], type='l', ylab = names(yieldDat)[i])
  print(paste(names(yieldDat)[i], "is", 
    if(adf.test(diffYield[, i])$p.value <= 0.01) "stationary" else "not stationary"))
}

pcaDiffYield = princomp(diffYield)
names(pcaDiffYield)
summary(pcaDiffYield)
plot(pcaDiffYield)

stocks = read.csv("F:\\R\\SaDAfFE\\data\\Stock_FX_Bond_2004_to_2006.csv")
ffData = read.table("F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 17\\FamaFrenchDaily.txt", header=T)

stockDates = as.POSIXlt(stocks[["DATE"]], '', '%m/%d/%Y')
stockDates = (stockDates$year + 1900)*10000 + (stockDates$mon + 1)*100 + stockDates$mday
ffDates = ffData[['date']]
commonDates = intersect(ffDates, stockDates)

stocks = stocks[stockDates %in% commonDates, ]
ffData = ffData[ffDates %in% commonDates, ]
ffData = ffData[-1,] # delete first row since stocks_diff lost a row due to differencing

stockSubset = as.matrix(stocks[, c("GM_AC","F_AC","UTX_AC","MRK_AC")])
stocksDiff = 100*apply(log(stockSubset), 2, diff) - ffData$RF
colnames(stocksDiff) = c("GM","Ford","UTX","Merck")
fit1 = lm(stocksDiff~ffData$Mkt.RF)
summary(fit1)

cor(fit1$resid)
cor.test(fit1$resid[, 1], fit1$residuals[, 2])
cor.test(fit1$resid[, 1], fit1$residuals[, 3])
cor.test(fit1$resid[, 2], fit1$residuals[, 3])

# variance estimation(Beta' %*% COV(Mkt.RF) %*% Beta) is very inaccurate.
fit1$coefficients[2,] %*% t(fit1$coefficients[2,]) * var(ffData$Mkt.RF)# + var(fit1$resid)
var(stocksDiff)

fit2 = lm(stocksDiff~Mkt.RF+SMB+HML, data=ffData)
summary(fit2)

cor(fit2$resid)
cor.test(fit2$resid[, 1], fit2$residuals[, 2])
cor.test(fit2$resid[, 1], fit2$residuals[, 3])
cor.test(fit2$resid[, 2], fit2$residuals[, 3])

# AIC does not support multiple response.
logLikGauss = function(vals){
  -length(vals) * 0.5 * (log(sum(vals^2)) + log(2*pi) + 1 - log(length(vals)))
}
myAIC = function(fit, method="aic"){
  coeficients = coef(fit)
  residual = residuals(fit)
  name = colnames(coeficients)
  aics = matrix(nrow=1, ncol=length(name), dimnames=list("1", name))
  if(method != "aic" & method != "bic") stop("only aic and bic are supported")
  isAIC = (method == "aic")
  for(i in 1:length(name)){
    p = length(coeficients[, name[i]])
    aics[1,i] = -2 * logLikGauss(residual[, i]) + (if(isAIC) 2 else log(length(residual[, i]))) * (p + 1)
  }
  if(isAIC) list(aic=aics) else list(bic=aics)
}
myAIC(fit1)
myAIC(fit2)
myAIC(fit1, "bic")
myAIC(fit2,"bic")

estVarByFactor = function(betas, factors, varResid){
  betas %*% cov(factors) %*% t(betas) +varResid
}

factorNames = c("Mkt.RF", "SMB", "HML")
stock1Beta = matrix(c(0.5, 0.4, -0.1), nrow=1, dimnames=list("1", factorNames))
stock1VarResid = 23
stock2Beta = matrix(c(0.6, 0.15, 0.7), nrow=1, dimnames=list("1", factorNames))
stock2VarResid = 37

factors = ffData[, factorNames]
estVarByFactor(stock1Beta, factors, stock1VarResid)
estVarByFactor(stock2Beta, factors, stock2VarResid)
# covriance between stock1 and stock2
stock1Beta %*% cov(factors) %*% t(stock2Beta)

stocks = read.csv("F:\\R\\SaDAfFE\\data\\Stock_FX_Bond.csv")
stocksAC = c("GM_AC", "F_AC", "UTX_AC", "CAT_AC", "MRK_AC", "PFE_AC", "IBM_AC", "MSFT_AC")
stocks = stocks[, stocksAC]
stocksLogReturn = apply(apply(stocks, 2, log), 2, diff)
factorAnalysis = factanal(stocksLogReturn, factors=2, rotation="none")
loads = matrix(as.numeric(loadings(factorAnalysis)), ncol=2)
rownames(loads) = stocksAC
(loads %*% t(loads))['F_AC', 'IBM_AC']
cor(stocksReturn[, 'F_AC'], stocksReturn[, 'IBM_AC'])

yields2009 = read.csv("F:\\R\\SaDAfFE\\data\\yields2009.csv",header=T)
pca = prcomp(yields2009[, -1])

data(midcapD.ts, package='fEcofin')
fact = factanal(midcapD.ts[,-1], factors=9)
summary(fact)
loadings(fact) %*% t(loadings(fact))
loadings(fact) %*% t(loadings(fact)) - cor(midcapD.ts[,-1])