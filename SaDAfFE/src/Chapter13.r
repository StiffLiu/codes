library(faraway)
set.seed(99)
len = 11
x = 1:len
x[len] = 50
y = 1 + x + rnorm(len)
y2 = c(y[1:(len - 1)], y[len] - 45)
x2 = c(x[1:(len - 1)], 5.5)

my_plot = function(x, y, l){
  main = paste("(", l, ")")
  plot(x, y, ylim=c(0, 60), cex=c(rep(1.25, 10), 1.5), pch=c(rep(21, 20), 19), main=main)
  model = lm(y~x)
  abline(model, lwd=2)
  plot(hatvalues(model), xlab="index", ylab="leverage", ylim=c(0,1), main=main)
  plot(rstudent(model), ylab="studentized leverage", main=main)
  plot(residuals(model), ylab="residual", main=main)
  plot(sqrt(cooks.distance(lm(y~x))), ylab="square root Cook's D", cex=1, main=main, ylim=c(0,11))
  halfnorm(sqrt(cooks.distance(lm(y~x))), ylab="square root Cook's D", cex=1, main=main, xlim=c(0,1.85))
}

par(mfrow=c(3, 6), lwd=1)
my_plot(x, y, "a")
my_plot(x, y2, "b")
my_plot(x2, y, "c")

dat = read.table(file="F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 13\\WeekIntAllData.txt",header=T)
attach(dat)
cm10_dif = diff(cm10)
aaa_dif = diff(aaa)
cm30_dif = diff(cm30)
ff_dif = diff(ff)
fit = lm(aaa_dif~cm10_dif+cm30_dif)
n=length(cm10)

par(mfrow=c(2,2))
plot(hatvalues(fit),ylab="Leverage",xlab="Index", main="(a)")
plot(2:n,rstudent(fit),ylab="rstudent",xlab="Index", main="(b)")
plot(2:n,cooks.distance(fit),ylab="Cook's D",xlab="Index", main="(c)")
plot(2:n,cooks.distance(fit),ylab="Cook's D", xlim=c(368,378),xlab="Index", main="(d)")

n = 80
set.seed("781235")
e = abs(matrix(runif(12*n),nrow=n) %*% rep(1,12))^4
e= e/mean(e) 
x1 = sort(runif(n))
x2 = rbeta(n,6,.5)
y =( 8*x2 + x1 + 5*x1^3) + ( 4* x2 + x1 + 7*x1^3) * e 

par(mfrow=c(1,2))
plot(x1,y,xlab=expression(x[1]))
plot(x2,y,xlab=expression(x[2]))

fit = lm(y~x1+x2)
rstudent = rstudent(fit)

par(mfrow=c(1,2))
qqnorm(rstudent,datax=T,main="Normal QQ Plot")
hist(rstudent,12)

par(mfrow=c(1,3))
plot(x1,rstudent,main="(a)",xlab=expression(x[1]))
fit2 = loess(rstudent~x1)
lines(x1,fit2$fitted)

plot(x2,rstudent,main="(b)",xlab=expression(x[1]))
fit3 = loess(rstudent~x2)
ordx2 = order(x2)
lines(x2[ordx2],fit3$fitted[ordx2])

fitquad = lm(y~poly(x1,2)+x2 )
rstudentquad = rstudent(fitquad)
plot(fitquad$fitted,abs(rstudentquad),xlab="fitted values",ylab="abs(rstudent)",main="(c) ")
fit4 = loess(abs(rstudentquad)~fitquad$fitted)
ord = order(fitquad$fitted)
lines(fitquad$fitted[ord],fit4$fitted[ord])

library(AER)
data(CPS1988)
attach(CPS1988)
# fitLm1 = lm(wage~education+experience+ethnicity)
fitLm1 = lm(log(wage) ~education+experience^2+ethnicity)

par(mfrow=c(3,2))
resid1 = rstudent(fitLm1)

plot(fitLm1$fit, resid1, main="(a)") # ylim=c(-1500, 1500), 
lines(lowess(fitLm1$fit, resid1, f=.2), lwd=5, col="red")
abline(h=0, col="blue", lwd=5)

plot(fitLm1$fit, abs(resid1), main="(b)") # ylim=c(0, 1500), 
lines(lowess(fitLm1$fit, abs(resid1), f=.2), lwd=5, col="red")
abline(h=mean(abs(resid1)), col="blue", lwd=5)

qqnorm(resid1, datax=F, main="(c)")
qqline(resid1, datax=F, lwd=5, col="blue")

plot(education, resid1, main="(d)") # ylim=c(-1000, 1500), 
lines(lowess(education, resid1), lwd=5, col="red")
abline(h=0, col="blue", lwd=5)

plot(experience, resid1, main="(e)") # ylim=c(-1000, 1500), 
lines(lowess(experience, resid1), lwd=5, col="red")
abline(h=0, col="blue", lwd=5)

library(faraway) # required for halfnorm
hatFit = hatvalues(fitLm1)
len = length(hatFit)
hatOrder = order(hatFit, (len:1), (1:len), decreasing=T)
sCooksDistFit = sqrt(cooks.distance(fitLm1))
len = length(sCooksDistFit)
cooksOrder = order(sCooksDistFit, (len:1), (1:len), decreasing=T)
len = length(wage)
wageOrder = order(wage, (len:1), (1:len), decreasing=T)
len = length(resid1)
residOrder = order(resid1, (len:1), (1:len), decreasing=T)
headCount = 3

my_plot = function(dat){
  plot(dat, ylab=deparse(substitute(dat)))
  for (i in 1:headCount){
    lines(c(hatOrder[i], hatOrder[i]), c(0, dat[hatOrder[i]]), col="red")
  }
  for (i in 1:headCount){
    lines(c(cooksOrder[i], cooksOrder[i]), c(0, dat[cooksOrder[i]]), col="blue")
  }
  for (i in 1:headCount){
    lines(c(wageOrder[i], wageOrder[i]), c(0, dat[wageOrder[i]]), col="yellow")
  }
  for (i in 1:headCount){
    lines(c(residOrder[i], residOrder[i]), c(0, dat[residOrder[i]]), col="green")
  }
}

par(mfrow=c(2,3))
my_plot(hatFit)
my_plot(sCooksDistFit)
my_plot(resid1)
#halfnorm(sCooksDistFit)
my_plot(wage)
my_plot(experience)
my_plot(education)