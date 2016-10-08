gasFlow = read.csv('F:\\R\\SaDAfFE\\data\\GasFlowData.csv')
f1 = gasFlow[[1]]
f2 = gasFlow[[2]]
useEmpericalCdf = FALSE
# Nonparamatric way: emperical cdf
npu1 = ecdf(f1)(f1)
npu2 = ecdf(f2)(f2)

library('fGarch')
# Paramatric way: skewed-t CDFs
# Fit using Maximum-Likelyhood-Estimator
sstdE1 = sstdFit(f1)$estimate
sstdE2 = sstdFit(f2)$estimate
sstu1 = psstd(f1, mean=sstdE1[1], sd=sstdE1[2], nu=sstdE1[3], xi=sstdE1[4])
sstu2 = psstd(f2, mean=sstdE2[1], sd=sstdE2[2], nu=sstdE2[3], xi=sstdE2[4])

if(useEmpericalCdf){
  u1 = npu1
  u2 = npu2
}else{
  u1 = sstu1
  u2 = sstu2
}
copEmprical = pempiricalCopula(npu1,npu2)

hist(u1)
hist(u2)
pairs(cbind(u1, u2))

z1 = qnorm(u1)
z2 = qnorm(u2)
pairs(cbind(z1, z2))

gasFlowMleEst = function(copula){
  result = fitCopula(copula, cbind(u1, u2), method="ml")
  contour(result@copula, pCopula,main= class(result@copula)[[1]])
  result
}

par(mfrow=c(2,3))
contour(copEmprical$x,copEmprical$y,copEmprical$z,main="Empirical")

# Kendall's tau correlation of the two flow data
# Use Kendall's tau correlation because it's invariant to strictly monnotonic transformations.
cor_tau = cor(u1,u2,method="kendall")

# The theory is that(for Gauss Copula(meta-Gauss distribution)):
#   KendallTau(X, Y) = 2 / pi * arcsin(Cor(X, y))
# So Cor(X,Y) = sin(pi / 2 * KendallTau(X, Y))
cor_start = sin(pi / 2 * cor_tau)
# Fit using Gauss Copula
gasFlowMleEst(normalCopula(cor_start))
# Fit using t Copula
gasFlowMleEst(tCopula(cor_start))
# Fit using Frank Copula, the starting parameter -1 is choosen randomly
gasFlowMleEst(frankCopula(-1))
# Fit using Clayton Copula, the starting parameter 2 is choosen randomly
gasFlowMleEst(claytonCopula(2))
# Not able to fit using Gumbel Copula, the data is kind of counter-monotonicity
#   but Gumbel Copula does not tends to be counter-monotonicity for any of it's parameter value.
# flowCCop = gasFlowMleEst(gumbelCopula(1.2))


numSeed = 5640
length = 1000
library(copula)
set.seed(numSeed)
cop_t_dim3 = tCopula(c(-.6,.75,0), dim=3, dispstr="un", df=3)
rand_t_cop = rCopula(length, cop_t_dim3)
cor(rand_t_cop)
cor.test(rand_t_cop[,1], rand_t_cop[, 2])
pairs(rand_t_cop)
# Pearson correlation changes with transformations
# cor(qt(rand_t_cop[,1], df=3), qt(rand_t_cop[,2], df=3)) ~= -.6

cop_normal_dim3 = normalCopula(c(-.6, .75, 0), dim = 3,  dispstr = "un")
mvdc_normal <- mvdc(cop_normal_dim3, c("exp", "exp", "exp"),
  list(list(rate=2), list(rate=3), list(rate=4)))
set.seed(numSeed)
rand_mvdc = rMvdc(length, mvdc_normal)
pairs(rand_mvdc)
par(mfrow=c(2,2))
plot(density(rand_mvdc[,1]))
plot(density(rand_mvdc[,2]))
plot(density(rand_mvdc[,3]))

library(Ecdat)  # need for the data
library(copula) # for copula functions
library(fGarch) # need for standardized t density
library(MASS) # need for fitdistr and kde2d
library(fCopulae) # additional copula functions (pempiricalCopula and ellipticalCopulaFit)
data(CRSPday, package="Ecdat")
ibm = CRSPday[,5]
crsp = CRSPday[,7]
est.ibm = as.numeric(fitdistr(ibm, "t")$estimate)
est.crsp = as.numeric(fitdistr(crsp, "t")$estimate)
est.ibm[2] = est.ibm[2] * sqrt(est.ibm[3] / (est.ibm[3] - 2))
est.crsp[2] = est.crsp[2] * sqrt(est.crsp[3] / (est.crsp[3] - 2))
cor_tau = cor(ibm, crsp, method="kendall")
# see line 47
omega = sin(pi / 2 * cor_tau)
cop_t_dim2 = tCopula(omega, dim=2, dispstr="un", df=4)

n = length(ibm)
# Paramatric way
data1 = cbind(pstd(ibm, mean=est.ibm[1], sd = est.ibm[2], nu=est.ibm[3]),
  pstd(crsp, mean=est.crsp[1], sd=est.crsp[2], nu=est.crsp[3]))
# Transform using rank
data2 = cbind(rank(ibm) / (n + 1), rank(crsp) / (n + 1))
ft1 = fitCopula(cop_t_dim2, method="mpl", optim.method="L-BFGS-B", data=data1,
  start=c(omega, 5), lower=c(0, 2.5), upper=c(.5, 15))
ft2 = fitCopula(cop_t_dim2, method="mpl", optim.method="L-BFGS-B", data=data2,
  start=c(omega, 5), lower=c(0, 2.5), upper=c(.5, 15))


mvdc_t_t = mvdc(cop_t_dim2, c("std", "std"),
  list(list(mean=est.ibm[1], sd=est.ibm[2], nu=est.ibm[3]),
       list(mean=est.crsp[1], sd=est.crsp[2], nu=est.crsp[3])))
start=c(est.ibm, est.crsp, ft1@estimate)
objFn = function(param){
  -loglikMvdc(param, cbind(ibm, crsp), mvdc_t_t)
}
t1 = proc.time()
fit_cop = optim(start, objFn, method="L-BFGS-B",
  lower = c(-.1, 0.001, 2.5, -.1, 0.001, 2.5, .2, 2.5),
  upper = c(.1, .03, 15, .1, .03, 15, .8, 15))
t2 = proc.time()
total_time = t2 - t1
total_time[3] / 60
mvdc_est = fit_cop$par
cop_est = tCopula(mvdc_est[7], df=as.integer(mvdc_est[8]))
mvdc_t_t_e = mvdc(cop_est, c("std", "std"),
  list(list(mean=mvdc_est[1], sd=mvdc_est[2], nu=mvdc_est[3]),
       list(mean=mvdc_est[4], sd=mvdc_est[5], nu=mvdc_est[6])))

myFitCopula = function(copulaFn, start=1){
  result = fitCopula(data=data1, copulaFn(start), optim.method="BFGS")
  contour(result@copula, pCopula,main= class(result@copula)[[1]])
  result
}
copEmprical = pempiricalCopula(data1[,1], data1[,2])
par(mfrow=c(3,2))

contour(copEmprical$x,copEmprical$y,copEmprical$z,main="Empirical")
contour(tCopula(param=ft2@estimate[1], df=as.integer(ft2@estimate[2])), pCopula, main="tCopula")
fnorm = myFitCopula(normalCopula, -.3)
fgumbel = myFitCopula(gumbelCopula, 1.1)
ffrank = myFitCopula(frankCopula)
fclayton = myFitCopula(claytonCopula)


  