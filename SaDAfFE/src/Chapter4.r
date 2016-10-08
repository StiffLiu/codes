logR = diff(log(EuStockMarkets))
plot(logR)
plot(as.data.frame(logR))

index.names = dimnames(logR)[[2]]
par(mfrow=c(2,2))
for(i in 1:4)
{
  qqnorm(logR[,i], datax=T, main=index.names[i])
  qqline(logR[,i], datax=T)
  print(shapiro.test(logR[,i]))
}

n = dim(logR)[1]
q.grid = (1:n)/(n + 1)
df = c(1, 4, 6, 10, 20, 30)
for(i in 1:4)
{
  windows()
  par(mfrow=c(3,2))
  for(j in 1:6)
  {
    qqplot(logR[,i], qt(q.grid, df=df[j]),
	  main=paste(index.names[i], ", df=", df[j]))
	abline(lm(qt(c(.25, .75), df=df[j])~quantile(logR[,i], c(.25, .75))))
  }
}

library("fGarch")
x = seq(-.1, .1, by=.001)
par(mfrow=c(1,1))
plot(density(logR[,1]), lwd=2, ylim=c(0, 8), xlim=c(0.02, 0.05))
lines(x, dstd(x, mean=median(logR[, 1]), sd=mad(logR[, 1]), nu=5), lty=5, lwd=2)
lines(x, dnorm(x, mean=mean(logR[, 1]), sd=sd(logR[,1])), lty=3, lwd=4)
legend("topleft", c("KDE", "t: df=5", "normal"), lwd=c(2, 2, 4), lty=c(1, 5, 3))

# Exercises 4.12.1
library("fEcofin")
fordReturn = ford.s[["FORD"]]
mean(fordReturn)
median(fordReturn)
sd(fordReturn)
plot(density(fordReturn))
qqnorm(fordReturn, datax=T)
qqline(fordReturn, datax=T)
print(shapiro.test(fordReturn))
df = c(1, 4, 5, 6, 8, 10, 20, 30)
par(mfrow=c(4,2))
for(i in 1:length(df))
{
  tData = rt(length(fordReturn), df=df[i])
  qqplot(fordReturn, tData, main=paste("Ford return df=", df[i]))
  abline(lm(quantile(tData, c(.25, .75))~quantile(fordReturn, c(.25, .75))))
}
sd(fordReturn)/sqrt(length(fordReturn))
prob = 0.5
qtl = quantile(fordReturn, prob)
dsty = approxfun(density(fordReturn))(qtl)
prob * (1 - prob)/(n * dsty * dsty)


# Exercises 4.12.2
library("Ecdat")
data()
dy = Garch[["dy"]]
plot(density(dy), lwd=2)
lines(density(rnorm(length(dy), mean=mean(dy), sd=mad(dy))), lty=5, lwd=2)
lines(density(rnorm(length(dy), mean=mean(dy), sd=sd(dy))), lty=3, lwd=4)
legend("topleft", c("KDE", "mean=mean, sd=mad", "mean=mean, sd=sd"), lwd=c(2, 2, 4), lty=c(1, 5, 3))


# Exercises 4.12.3
library("Ecdat")
data()
diffbp = diff(Garch[["bp"]])
diffbp = diff(log(Garch[["bp"]]))
pl = c(.25, .1, .05, .025, .01, .0025)
par(mfrow=c(3,2))
for(i in 1:length(pl))
{
  nData = rnorm(length(diffbp))
  p = pl[i]
  qqplot(diffbp, nData, main=paste("Dollar to pound, ref line p = ", pl[i]))
  abline(lm(quantile(nData, c(p, 1 - p))~quantile(diffbp, c(p, 1 - p))))
}