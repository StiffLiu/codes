dat = read.csv("F:\\R\\SaDAfFE\\data\\Stock_FX_Bond.csv", header=T)
prices = cbind(dat$GM_AC, dat$F_AC, dat$CAT_AC, dat$UTX_AC, dat$MRK_AC, dat$IBM_AC)
n = dim(prices)[1]
d = dim(prices)[2]
returns = 100 * (prices[2:n,] / prices[1:(n-1), ] - 1)
pairs(returns)
mean_vect = apply(returns, 2, mean)
cov_mat = cov(returns)
sd_vect = sqrt(diag(cov_mat))
library(quadprog)
Amat = cbind(rep(1, d), mean_vect, diag(1, d), diag(-1, d)) #equality constratis are at the beginning
maxR = 0.0815
muP = seq(.04,maxR,length=300)
weights = matrix(ncol=d)
sds = c()
for(i in 1:length(muP)){
  bvec = c(1, muP[i], rep(-0.1, d), rep(-0.5, d))
  result = solve.QP(Dmat = 2 * cov_mat, dvec = rep(0, d), Amat=Amat, bvec=bvec, meq = 2)
  sds = c(sds, result$value)
  weights = rbind(weights, result$solution)
}
rfRate = 3.0 / 365
sharpe = (muP - rfRate) / sds
tangencyIndex = (sharpe == max(sharpe))
tangency = sharpe[tangencyIndex]
efficientIndex = (muP > muP[tangencyIndex])
plot(sds, muP, type='l')
lines(c(0, (maxR - rfRate) / tangency), c(rfRate, maxR))
lines(sds[efficientIndex], muP[efficientIndex], type='l', col="red")
lines(c(0, sds[tangencyIndex]), c(rfRate, muP[tangencyIndex]), type='l', col="red")
points(c(sds[tangencyIndex]), c(muP[tangencyIndex]))

# Is there any way to find the index?
index0_07 = (1:length(muP))[(abs(muP - 0.07) == min(abs(muP - 0.07)))]
weights[index0_07, ]