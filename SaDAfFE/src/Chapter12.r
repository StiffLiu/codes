dat = read.table(file="F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 12\\WeekInt.txt",header=T)
attach(dat)
aaa_dif = diff(aaa)
cm10_dif = diff(cm10)
ff_dif = diff(ff)
cm30_dif = diff(cm30)
fit1 = lm(aaa_dif~cm10_dif + cm30_dif + ff_dif)
anova(fit1)
anova(lm(aaa_dif~cm10_dif), fit1)
anova(lm(aaa_dif~cm10_dif + cm30_dif), fit1)

subsets = regsubsets(aaa_dif~.,
  data=as.data.frame(cbind(cm10_dif,cm30_dif,ff_dif)),nbest=1)
sumSubset = summary(subsets)                      #  Example 12.7
sumSubset

par(mfrow=c(1,3),lab=c(2,5,3),pch=19)
plot(1:3,sumSubset$bic,type="b",xlab="number of variables",
   ylab="BIC",cex=2.5)
plot(1:3,sumSubset$cp,type="b",xlab="number of variables",
   ylab="Cp",cex=2.5)
plot(1:3,sumSubset$adjr2,type="b",xlab="number of variables",
   ylab="adjusted R2",cex=2.5)

library(faraway)
vif(lm(aaa_dif~cm10_dif+cm30_dif+ff_dif))   #  Example 12.8

data(nelsonplosser, package='fEcofin')
names(nelsonplosser)
new_np = na.omit(nelsonplosser)
n = dim(new_np)[1]
attach(new_np)
year = X.Y.m.d[2:length(X.Y.m.d)]
par(mfrow=c(2,3),cex.lab=1.35)
plot(year,diff(gnp.r),type="b",ylab="differences",main="gnp.r")
plot(year,diff(log(gnp.r)),type="b",ylab="differences",main="log(gnp.r)")
plot(year,diff(sqrt(gnp.r)),type="b",ylab="differences",main="sqrt(gnp.r)")
plot(year,diff(ip),type="b",ylab="differences",main="ip")
plot(year,diff(log(ip)),type="b",ylab="differences",main="log(ip)")
plot(year,diff(sqrt(ip)),type="b",ylab="differences",main="sqrt(ip)")
lmReg = lm(formula = diff(log(sp)) ~ diff(gnp.r) + diff(gnp.pc) + diff(log(ip)) + diff(log(cpi)) + diff(emp) + diff(bnd), data = new_np)
summary(lmReg)
anova(lmReg)
vif(lmReg)
stepAIC(lm(formula = diff(log(sp)) ~ diff(gnp.r) + diff(gnp.pc) + diff(log(ip)) + diff(log(cpi)) + diff(emp) + diff(bnd)))

library(leaps)
predictors = cbind(gnp.r, gnp.pc, log(ip), log(cpi))
colnames(predictors) = c('gnp.r', 'gnp.pc', 'log(ip)', 'log(cpi)')
cpResults = leaps(cbind(gnp.r, gnp.pc, log(ip), log(cpi)), log(sp))
selections = cpResults$which
selections[selections] = 1
selections = cbind(selections, cpResults$Cp)
colnames(selections) = c('gnp.r', 'gnp.pc', 'log(ip)', 'log(cpi)', 'cp')

fit1 = lm(aaa_dif ~ cm10_dif)
fit2 = lm(aaa_dif~cm10_dif+cm30_dif)
fit4 = lm(aaa_dif~cm30_dif)

par(mfrow=c(2,2))
crPlot(fit2,var="cm10_dif",main="(a)",smooth=F,lty=1,lwd=2,col="black")
crPlot(fit2,var="cm30_dif",main="(b)",smooth=F,lty=1,lwd=2,col="black")
plot(cm10_dif,aaa_dif,main="(c)")
regLine(fit1,col="black",lwd=2,lty="dashed")
plot(cm30_dif,aaa_dif,main="(d)")
regLine(fit4,col="black",lwd=2,lty="dashed")

library(MASS)
library(AER)
data(USMacroG)
attach(as.data.frame(USMacroG))
MacroDiff = apply(USMacroG, 2, diff)
pairs(cbind(consumption, dpi, cpi, government, unemp))
fitLm1 = lm(consumption~dpi+cpi+government+unemp)
summary(fitLm1)
confint(fitLm1)
fitLm2 = stepAIC(fitLm1)
summary(fitLm2)
AIC(fitLm1)-AIC(fitLm2)