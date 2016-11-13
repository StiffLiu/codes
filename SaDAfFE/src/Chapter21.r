library(mgcv)
data(CPS1988, package="AER")
attach(CPS1988)
fitGam = gam(log(wage)~s(education)+s(experience)+ethnicity)
summary(fitGam)
par(mfrow=c(1,2))
plot(fitGam)

# The CKLS model for interest rate: µ(t, r) = a(t) {θ(t) − r}.
# Where a is a liear function of time "t"
#  θ(t) is a piecewise linear spline with knots.
# "r" is the lagged rate
data(Irates, package="Ecdat")
r1 = Irates[,1]
n = length(r1)
lagr1 = lag(r1)[-n]
diffr1 = diff(r1)
n = length(lagr1)
knots = seq(from=1950, to=1985, length=10)
t = seq(from=1946, to=1991+2/12, length=n)
X1 = outer(t, knots, FUN="-")
X2 = X1 * (X1 > 0)

# X3 is a matrix, whose colum is a vector:
#   1. the plus function: (t - Ki)+ (for i > 1)
#      Ki is 1946 (for i = 2)
#      Ki is the knots vector (for i > 2)
#   2. all ones (for i = 1)
X3 = cbind(rep(1, n), (t - 1946), X2)
m2 = dim(X3)[2]
m = m2 - 1

# a_t = X3 %*% a
# theta_t = X3 %*% theta
nlmod_CKLS_ext = nls(diffr1~X3[,1:2] %*% a * (X3 %*% theta - lagr1),
  start=list(theta=c(10, rep(0, m)), a=c(0.01, 0)),
  control=list(maxiter=200))
AIC(nlmod_CKLS_ext)
param4 = summary(nlmod_CKLS_ext)$parameters[,1]
theta = param4[1:m2]
a = param4[(m2+1):length(param4)]
a_t = X3[,1:2] %*% a
theta_t = X3 %*% theta
res = residuals(nlmod_CKLS_ext)
resSquare = res^2
nlmod_CKLS_ext_res = nls(resSquare ~ A * lagr1^B, start=list(A=.2, B=0.5))

par(mfrow=c(1,2))
plot(t, theta_t, ylim=c(0,16), ylab="rate", main="(a)",
  col="red", type="l", lwd=2)
lines(t, lagr1)
lines(t, a_t, col="blue", lwd=1)
lines(t, a_t * (theta_t - lagr1) + lagr1, col="green", lwd=1)
legend("topleft", c("theta(t)", "a(t)", "a(t)*(theta(t) - lagr) + lagr", 
  "lagged rate"), lwd=c(2, 1, 1, 1), col=c("red", "blue", "green", "black"))
plot(lagr1, abs(res), pch=5, ylim=c(0, 6), ylab="", main="(b)")
lines(lagr1, sqrt(fitted(nlmod_CKLS_ext_res)), lw=3, col="red", type="l")
legend("topleft", c("abs res", "volatility fn"),
  lty=c(NA, 1), pch=c(5, NA), col=c("black", "red"), lwd=c(1, 2))