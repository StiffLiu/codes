data(Icecream, package='Ecdat')
attach(Icecream)
fit_ic_lm = lm(formula=cons~income + price + temp, data=Icecream)
library(car)
durbinWatsonTest(fit_ic_lm)
fit_ic_ar_lm = arima(x=cons, order=c(1, 0, 0), xreg=cbind(income, price, temp))
fit_ic_ma_lm = arima(x=cons, order=c(0, 0, 1), xreg=cbind(income, price, temp))
summary(fit_ic_ar_lm)
summary(fit_ic_ma_lm)
par(mfrow=c(1, 3))
acf(fit_ic_lm$residuals, main="linear model/indep.noise")
acf(fit_ic_ar_lm$residuals, main="linear model/AR(1) noise")
acf(fit_ic_ma_lm$residuals, main="linear model/MA(1) noise")

x = -10:10
X = cbind(rep(1, length(x)), x)
B = X %*% solve(t(X) %*% X)
Phi = seq(-0.75, 0.75, 0.25)
covBeta = matrix(nrow=length(Phi), ncol=2)
for(i in 1:length(Phi)){
  covError = toeplitz(Phi[i]^(0:(length(x) - 1)))
  covBeta[i, ] = sqrt(diag(t(B) %*% covError %*% B))
}

plot(Phi, covBeta[4,1]/covBeta[, 1], type="b", 
  pch="o", xlim=c(-.75, .75), ylim=c(0.4, 2.5), ylab="SE ratio")
par(new=T)
plot(Phi, covBeta[4,2]/covBeta[, 2], type="b", ylab="",
  pch="*", lty=3, xlim=c(-.75, .75), ylim=c(0.4, 2.5))
legend("topright", c("intercept", "slop"), 
  pch=c("o", "*"), lty=c(1, 3))

bondprices = read.table("F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 14\\bondprices.txt",header=T)
fit = nls(price~1000 * exp(-r * maturity), data=bondprices, start=list(r=.04))
summary(fit)

attach(bondprices)
plot(maturity, price - summary(fit)$resid, type="l")
points(maturity, price, pch="*")
legend("topright", c("price", "predicted price"), pch=c("*", NA), lty=c(NA, 1))

defaultFreq = read.table("F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 14\\DefaultData.txt",header=T)
attach(defaultFreq)
freq = freq/100

fit_nls = nls(freq~exp(b0 + b1 * rating), start=list(b0=1, b1=1))
resid_nls = summary(fit_nls)$resid
fitted_nls = freq - resid_nls

fit_bow = lm(log(freq[freq > 0])~rating[freq > 0])
resid_bow = summary(fit_bow)$resid
fitted_bow = log(freq[freq > 0]) - resid_bow

fit_tbs = nls(sqrt(freq)~exp((b0 + b1 * rating)/2), start=list(b0=1, b1=1))
resid_tbs = summary(fit_tbs)$resid
fitted_tbs = sqrt(freq) - resid_tbs
tbs_loess = loess(abs(resid_tbs)~fitted_tbs, span=1, deg=1)

par(mfrow=c(1,2))
plot(fitted_tbs, abs(resid_tbs), type="p")
lines(fitted_tbs, tbs_loess$fit)
qqnorm(resid_tbs, datax=T, main="")
qqline(resid_tbs, datax=T)

stripPrices = read.table("F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 14\\strips_dec95.txt",header=T)
attach(stripPrices)
orderT = order(T)
T = T[orderT]
price = price[orderT]
forwardEmpirical = - diff(log(price)) / diff(T)
par(mfrow=c(1, 2))
plot(T, price, xlab="maturity", main="(a)")
plot(T[2:length(T)], forwardEmpirical, xlab="maturity", ylab="empirical forwrad-rate", type="b", main="(b)")

dataMatrix = rep(1, length(T))

# forwardRateCurve = forwardRate ~ beta0*T + (beta1*T^2)/2
# priceCurve = price~100*exp(-forwardRateCurve)
dataMatrix = cbind(dataMatrix, T)
priceCurve = price~100*exp(-(beta0*T + (beta1*T^2)/2))
fitLin = nls(priceCurve, data=stripPrices, start=list(beta0=0.03, beta1=0))
sumLin = summary(fitQuad)
forwardLin = dataMatrix %*% sumLin$coef[,1]

# forwardRateCurve = forwardRate ~ beta0*T + (beta1*T^2)/2 + (beta2*T^3)/3
# priceCurve = price~100*exp(-forwardRateCurve)
dataMatrix = cbind(dataMatrix, T^2)
priceCurve = price~100*exp(-(beta0*T + (beta1*T^2)/2 + (beta2*T^3)/3))
fitQuad = nls(priceCurve, data=stripPrices, start=list(beta0=0.03, beta1=0, beta2=0))
sumQuad = summary(fitQuad)
forwardQuad = dataMatrix %*% sumQuad$coef[,1]

# forwardRateCurve = forwardRate ~ beta0*T + (beta1*T^2)/2 + (beta2*T^3)/3 + (T>15)*(beta3*(T-15)^3)/3
# priceCurve = price~100*exp(-forwardRateCurve)
dataMatrix = cbind(dataMatrix, (T>15)*(T-15)^2)
priceCurve = price~100*exp(-(beta0*T + (beta1*T^2)/2 + (beta2*T^3)/3 + (T>15)*(beta3*(T-15)^3)/3))
fitSpline = nls(priceCurve, data=stripPrices, start=list(beta0=0.03, beta1=0, beta2=0, beta3=0))
sumSpline = summary(fitSpline)
forwardSpline = dataMatrix %*% sumSpline$coef[,1]

# forwardRateCurve = forwardRate ~ beta0*T + (beta1*T^2)/2 + (beta2*T^3)/3 + (beta3*T^4)/4
# priceCurve = price~100*exp(-forwardRateCurve)
dataMatrix[,4] = T^3
priceCurve = price~100*exp(-(beta0*T + (beta1*T^2)/2 + (beta2*T^3)/3 + (beta3*T^4)/4))
fitCubic = nls(priceCurve, data=stripPrices, start=list(beta0=0.03, beta1=0, beta2=0, beta3=0))
sumCubic = summary(fitCubic)
forwardCubic = dataMatrix %*% sumCubic$coef[,1]

plot(T[2:length(T)], forwardEmpirical, type="p", pch="*")
lines(T, forwardQuad, lty=1)
lines(T, forwardCubic,lty=3, lwd=3)
lines(T, forwardSpline,lty=2, lwd=2)
legend("bottomleft", c("quaratic", "cubic", "spline", "empirical"),
  pch=c(NA, NA, NA, "*"), lty=c(1, 3, 2, NA), lwd=c(1, 3, 2, NA))
  
####################
n = 80
set.seed("781235")
e = matrix(runif(12*n),nrow=n) %*% rep(1,12)
e = abs(e)^4
e= e/mean(e) 
x1 = runif(n)
x1 = sort(x1) 
x2 = rbeta(n,6,.5)
y =( 8*x2 + x1 + 5*x1^3) + ( 4* x2 + x1 + 7*x1^3) * e
library(MASS)
boxcox(lm(y~x1 + x1^2 + x2))
####################

par(mfrow=c(3,2))
plot(rating, freq - resid_nls, type="l", ylab="frequency")
points(rating, freq, pch="*")
legend("topleft", c("exponential", "data"), pch=c(NA, "*"), lty=c(1, NA))
summary(fit_nls)

plot(rating, log(freq + 1e-6), pch="o", ylab="log(default probability)")
lines(rating[freq > 0], fitted_bow)
lines(rating, log(fitted_nls), lty=2)
lines(rating, 2 * log(fitted_tbs), lty=3)
legend("topleft", c("bow", "nonlinear", "tbs", "data"),
  pch=c(NA, NA, NA, "o"), lty=c(1, 2, 3, NA))

fitted_resid = function(fitted_value, resid_value){ 
  fitted_value_o = order(fitted_value)
  fit_loess  = loess(abs(resid_value)~ fitted_value,span=1,deg=1)
  plot(fitted_value[fitted_value_o], abs(resid_value[fitted_value_o]), ylab="absolute residual")
  lines(fitted_value[fitted_value_o], fit_loess$fit[fitted_value_o])
  qqnorm(resid_value, datax=T)
  qqline(resid_value, datax=T)
}

fitted_resid(fitted_nls, resid_nls)
fitted_resid(fitted_tbs, resid_tbs)

library(AER)
data(CreditCard)
names(CreditCard)
attach(CreditCard)
validIndex = (age > 1)
CreditCard_clean = data.frame(card=card[validIndex], reports=reports[validIndex],
  income=income[validIndex], share=share[validIndex], age=age[validIndex],
  owner=owner[validIndex], dependents=dependents[validIndex], months=months[validIndex],
  majorcards=majorcards[validIndex])
attach(CreditCard_clean)
par(mfrow=c(3,3))
hist(dependents, main="dependents")
hist(income, main="income")
hist(as.integer(owner), main="owner")
hist(age, main="age")
hist(share, main="share")
hist(reports, main="reports")
hist(months, main="months")
hist(log(share), main="log(share)")
hist(log(reports + 1), main="log(reports + 1)")

fit1 = glm(formula = card ~ log(reports + 1) + income + 
  log(share) + age + owner + dependents + months,
  family = "binomial", data = CreditCard_clean)

library(MASS)
bestFit = stepAIC(fit1)
summary(bestFit)

# use the corrected data: attach(CreditCard_clean)
center = function(x){x - mean(x)}
log_reports_c = center(log(reports + 1))
income_c = center(income)
log_share_c = center(log(share))
dependents_c = center(dependents)

centeredFit = glm(formula=card~log_reports_c + income_c + log_share_c
  + dependents_c, family = "binomial")
summary(centeredFit)
# The fitted value are probabilities:
#  1/(1 + exp(X %*% beta))
 1/(1 + exp(-cbind(rep(1,length(log_reports_c)), log_reports_c, income_c, log_share_c, dependents_c) %*% centeredFit$coefficients)) - centeredFit$fitted


my_plot = function(x, component, ...){
  len = length(x)
  # centeredFit$centeredFit$coefficients:
  X = matrix(nrow=len, ncol=length(centeredFit$coefficients))
  #    Beta0 + Beta1 * reports + Beta2 * income + Beta3 * share + Beta4 * dependents
  indices = list(reports=2, income=3, share=4, dependents=5)
  index = indices[[component]]
  if(is.null(index)) stop(paste("unsupported ", component))
  l_const = rep(1, len)
  l_reports = if(index == 2) log(x + 1) - mean(log(reports + 1)) else rep(mean(log_reports_c), len)
  l_income = if(index == 3) x - mean(income) else rep(mean(income_c), len)
  l_share = if(index == 4) x - mean(log(share)) else rep(mean(log_share_c), len)
  l_dependents =if(index == 5) x - mean(dependents) else rep(mean(dependents_c), len)
  l_X = cbind(l_const, l_reports, l_income, l_share, l_dependents)
  l_Beta = centeredFit$coefficients
  y = 1/(1 + exp(-l_X %*% l_Beta))
  plot(x, y, ...)
}
par(mfrow=c(2,2))
my_plot(0:14, "reports", ylim=c(0, 1), xlab="reports", ylab="P(accept)")
my_plot(0:14, "income", ylim=c(0, 1), xlab="income", ylab="P(accept)")
my_plot(seq(-9, 0, 0.1), "share", ylim=c(0, 1), xlab="log(share)", ylab="P(accept)")
my_plot(0:6, "dependents", ylim=c(0, 1), xlab="dependents", ylab="P(accept)")
par(mfrow=c(2,2))
plot(log(share), card)
plot(log(share), reports)
plot(log(share), income)
plot(log(share), majorcards)

beta1 = 1
beta2 = -1
sigma = 0.02
x = seq(-1, 2.5, 0.025)
y = beta1 * exp(beta2 * x) + rnorm(length(x), sd = sigma)
fit1 = lm(formula=log(y)~x)
par(mfrow=c(2,2))
plot(x,y)
plot(x, log(y))
qqnorm(fit1$residuals)
qqline(fit1$residuals)
plot(fit1$fitted, fit1$residuals)

x=seq(0,10,by=.4)
n = length(x)
y = 2 + 5*x + rnorm(n)
ind = c(3,24)
y[ind[1]] = y[ind[1]] + 35
y[ind[2]] = y[ind[2]] - 15

library(robust)
fit1 = lm(y~x)
robFit = ltsReg(y~x)
plot(x[-ind], y[-ind])
points(x[ind], y[ind], pch='*')
abline(fit1)
abline(robFit, lty=5)
legend("bottomright", c("good data", "residuals outlier", 
  "with outlier", "withut outlier"), 
  pch=c('o', '*', NA, NA), lty=c(NA, NA, 1, 5))

set.seed(99)
len = 11
x = 1:len
x[len] = 50
y = 1 + x + rnorm(len)
y2 = c(y[1:(len - 1)], y[len] - 45)
x2 = c(x[1:(len - 1)], 5.5)
my_plot = function(x, y, ...){
  plot(x, y, pch='*', cex=2, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), ...)
  abline(lm(y~x), col="red", lty=3)
  abline(ltsReg(y~x))
}
par(mfrow=c(1,3))
my_plot(x, y)
my_plot(x, y2, ylim=c(0, 50))
my_plot(x2, y)

library(AER)
data(USMacroG)
MacroDiff = as.data.frame(apply(USMacroG, 2, diff))
attach(MacroDiff)
fitNAR1 = arima(unemp, order=c(1, 0, 0), xreg=cbind(invest, government))
fitNAR2 = arima(unemp, order=c(2, 0, 0), xreg=cbind(invest, government))
fitNAR1MA1 = arima(unemp, order=c(1, 0, 1), xreg=cbind(invest, government))
fitLM = lm(formula=unemp~invest+government)
AIC(fitNAR1)
AIC(fitNAR2)
AIC(fitNAR1MA1)
AIC(fitLM)
BIC(fitNAR1)
BIC(fitNAR2)
BIC(fitNAR1MA1)
BIC(fitLM)

par(mfrow=c(2,2))
acf(fitNAR1$resid, main="arima(1, 0, 0)")
acf(fitNAR2$resid, main="arima(2, 0, 0)")
acf(fitNAR1MA1$resid, main="arima(1, 0, 1)")
acf(fitLM$resid, main="linear regression")

par(mfrow=c(2,2))
qqnorm(fitNAR1$resid, main="arima(1, 0, 0)")
qqline(fitNAR1$resid)
qqnorm(fitNAR2$resid, main="arima(2, 0, 0)")
qqline(fitNAR2$resid)
qqnorm(fitNAR1MA1$resid, main="arima(1, 0, 1)")
qqline(fitNAR1MA1$resid)
qqnorm(fitLM$resid, main="linear regression")
qqline(fitLM$resid)

library(Ecdat)
data(Irates)
r1 = Irates[, 1]
lag_r1 = lag(r1)[-length(r1)]
delta_r1 = diff(r1)
n = length(lag_r1)
par(mfrow=c(3, 2))
plot(r1, main="(a)")
plot(delta_r1, main="(b)")
plot(delta_r1^2, main="(c)")
plot(lag_r1, delta_r1, main="(d)")
plot(lag_r1, delta_r1^2, main="(e)")

nlmod_CKLS = nls(delta_r1 ~ a * (theta - lag_r1),
  start = list(theta=5, a=0.01), control = list(maxiter=200))
param = summary(nlmod_CKLS)$parameters[,1]
par(mfrow=c(2,2))
plot(lag_r1, ylim=c(0, 16), ylab="rate and theta",
  main="(a)", type="l")
abline(h=param[[1]], lwd=2, col="red")

res_sq = residuals(nlmod_CKLS)^2
nlmod_CKLS_res = nls(res_sq ~ A*lag_r1^B, start=list(A=.2, B=.5))
param2 = summary(nlmod_CKLS_res)$parameters[,1]
plot(lag_r1, sqrt(res_sq), pch=5, ylim=c(0,6), main="(b)")
attach(as.list(param2))
curve(sqrt(A*x^B), add=T, col="red", lwd=3)

nlmod_CKLS_wt = nls(delta_r1 ~ a * (theta - lag_r1),
  start=list(theta=5, a=.01), control=list(maxiter=200),
  weights=1/fitted(nlmod_CKLS_res))
plot(lag_r1, ylim=c(0, 16), ylab="rate and theta",
  main="(c)", type="l")
param3 = summary(nlmod_CKLS_wt)$parameters[,1]
abline(h=param3[1], lwd=2, col="red")

library(AER)
data(HousePrices)
fit1 = lm(price~., data=HousePrices)
summary(fit1)

library(MASS)
fit2 = boxcox(fit1, xlab=expression(alpha))

alpha = fit2$x[(1:length(fit2$y))[fit2$y==max(fit2$y)]]
transPrice = box.cox(HousePrices[["price"]], alpha)

library(car)
THousePrices = HousePrices
THousePrices[["price"]] = NULL
fit3 = lm(transPrice~.,data=THousePrices)
summary(fit3)
AIC(fit3)
BIC(fit3)

fit1 = glm(aircon~., family="binomial", data=HousePrices)
summary(fit1)
fit2 = stepAIC(fit1)
summary(fit2)

# Nelson-Siegel Forward rate
nsfr = function(theta0, theta1, theta2, theta3, tm){
  theta0 + (theta1 + theta2 * tm) * exp(-theta3 * tm)
}

# Nelson-Siegel Yield to maturity
nsy2m = function(theta0, theta1, theta2, theta3, tm){
  t2_dv_by_t3 = theta2 /theta3
  exp_neg_t3_mul_tm = exp(-theta3 * tm)
  theta0*tm + (theta1 + t2_dv_by_t3) * (1 - exp_neg_t3_mul_tm) / theta3 - t2_dv_by_t3 * tm * exp_neg_t3_mul_tm   
}

# Nelson-Siegel Yield
nsy = function(theta0, theta1, theta2, theta3, tm){
  t2_dv_by_t3 = theta2/theta3
  t3_mul_tm = theta3 * tm
  exp_neg_t3_mul_tm = exp(-t3_mul_tm)
  theta0 + (theta1 + t2_dv_by_t3) * (1 - exp_neg_t3_mul_tm) / t3_mul_tm
    - t2_dv_by_t3 * exp_neg_t3_mul_tm   
}

# Discount for Nelson-Siegel Yield
dnsy2m = function(theta0, theta1, theta2, theta3, tm){
  exp(-nsy2m(theta0, theta1, theta2, theta3, tm))
}

# Nelson-Siegel, price residuals
nsr = function(theta0, theta1, theta2, theta3, tm, price, par){
  price - par * dnsy2m(theta0, theta1, theta2, theta3, tm)
}

# log likelihood for i.i.d gauss distribution
logLikGauss = function(vals){
  -length(vals) *(log(sd(vals) * sqrt(2*pi)) + 0.5)
}

# Nelson-Siegel model log likelihood assuming residuals are i.i.d gauss variables
nsllgr = function(theta0, theta1, theta2, theta3, tm, price, par){
  logLikGauss(nsr(theta0, theta1, theta2, theta3, tm, price, par))  
}

# Nelson-Siegel model, sum of squares for residuals
nssosr = function(theta0, theta1, theta2, theta3, tm, price, par){
  sum(nsr(theta0, theta1, theta2, theta3, tm, price, par) ^ 2)
}

# Guess the initial value for Nelson-Siegel model
# This is done by assuming theta3 is one and lm to fit accordingly
# Then calculate parameter values based on lm fit result
ns_init_guess = function(price, maturity, par_value){
  minus_log_discount = -log(price/par_value)
  exp_neg_maturity = exp(-maturity)
  comp0 = 1 - exp_neg_maturity
  comp1 = maturity * exp_neg_maturity
  # theta3 = 1: Intecept, theta0, theta1 + theta2, -theta2
  tmp_theta = coef(lm(minus_log_discount~maturity+comp0+comp1))
  c(tmp_theta[2], tmp_theta[3] - tmp_theta[4], -tmp_theta[4], 1)
}

nsp = function(theta0, theta1, theta2, theta3, tm, par){
  par * dnsy2m(theta0, theta1, theta2, theta3, tm)
}

dat = read.table("F:\\R\\SaDAfFE\\data\\ZeroPrices.txt",header=T)
attach(dat)
par_value = 100
my_func = function(theta, func){
  func(theta[1], theta[2], theta[3], theta[4], maturity, price, par_value)
}

my_dnsy2m = function(param){
 if(!is.null(param$par)) param = param$par
 dnsy2m(param[1], param[2], param[3], param[4], maturity)
}

my_nsp = function(param){
 if(!is.null(param$par)) param = param$par
 nsp(param[1], param[2], param[3], param[4], maturity, par_value)
}

my_nsllgr = function(theta){-my_func(theta, nsllgr)}
my_nssosr = function(theta){my_func(theta, nssosr)}

my_inital_guess = function(){ ns_init_guess(price, maturity, par_value)}

theta_start = my_inital_guess()
est_sosr = optim(theta_start, my_nssosr)
est_llgr = optim(theta_start, my_nsllgr)
my_nssosr(est_sosr$par)
my_nssosr(est_llgr$par)
my_nsllgr(est_sosr$par)
my_nsllgr(est_llgr$par)
llgr_fit = my_nsp(est_llgr)
sosr_fit = my_nsp(est_sosr)
# This fit was produced with a bug in the code when writing.
another_fit = nsp(0.04292943,  0.27639571, -0.18983157,  1.05293383, maturity, par_value)
llgr_resid = price - llgr_fit
sosr_resid = price - sosr_fit
another_resid = price - another_fit

par(mfrow=c(2,2))
plot(price, type="p", pch="*")
lines(llgr_fit, col="red")
lines(sosr_fit, col="blue")
lines(another_fit, col="green")
legend("topright", c("sample", "log likelihood", "sum of squares", "another fit"), 
  lty=c(NA, 1, 1, 1), pch=c('*', NA, NA, NA), col=c("black", "red", "blue", "green"))
my_qq = function(dat, desc){
  qqnorm(dat, main=desc)
  qqline(dat)
}
my_qq(llgr_resid, "log likelihood")
my_qq(sosr_resid, "sum of squares")
my_qq(another_fit, "another fit")