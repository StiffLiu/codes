dat = read.csv("F:\\R\\SaDAfFE\\data\\Stock_FX_Bond_2004_to_2006.csv", header=T)
marketAC = "SP_AC"
riskFreeCol = "treasury"
stockAC = c("GM_AC", "F_AC", "UTX_AC", "CAT_AC", "MRK_AC", "PFE_AC", "IBM_AC")
cols = c(stockAC, marketAC)
prices = dat[, cols]
n = dim(prices)[1]

# Convert to daily percentage returns.
dat2 = as.matrix(cbind(dat[(2:n),"Three_month_treasury"]/365, 100*(prices[2:n,]/prices[1:(n-1),] - 1)))
colnames(dat2)[1] = riskFreeCol
riskFree = dat2[,riskFreeCol]
extRet = dat2[,2:ncol(dat2)] - riskFree
marketExtRet = extRet[, marketAC]
stockExtRet = extRet[, stockAC]

fitReg = lm(stockExtRet~marketExtRet)
fitReg1 = lm(stockExtRet~marketExtRet-1)
summary(fitReg)
res = residuals(fitReg)
pairs(res)
options(digits=6)
betas = fitReg$coeff[2,]
betas1 = fitReg1$coeff[1,]

apply(stockExtRet,2,mean)
betas*mean(marketExtRet)

fittedVar = diag(var(res)) + (betas^2)*var(marketExtRet)
# By CAPM, the residuals regarding the Market Portfolio of two different stocks
# should be uncorrelated, so ignnore the correlation of two different stocks
# the value on the off diagonal elements of fittedCor should be similiar to cor(stockExtRet)
marketCov = (betas %*% t(betas)) * var(marketExtRet)
varianceProduct = sqrt(fittedVar %*% t(fittedVar))
fittedCor = marketCov / varianceProduct
#variance due to market 
diag(marketCov) / fittedVar
# actualCor should be the same as cor(stockExtRet)
actualCor = fittedCor +  + var(res) / varianceProduct