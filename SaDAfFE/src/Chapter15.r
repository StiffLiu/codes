library(tseries)
library(urca)

yieldDat = read.table("F:\\R\\SaDAfFE\\data\\RPrograms\\Chapter 15\\treasury_yields.txt", header=T)
dat = as.matrix(yieldDat[, 3:7])
year = as.Date(yieldDat[["Date"]], "%m/%d/%y")
residual = residuals(lm(dat[,3]~dat[,1]+dat[,2]+dat[,4]+dat[,5]))

par(mfrow=c(2, 4))
mains = c("3-month", "6-month", "1-year", "2-year", "3-year")
for(i in 1:length(mains)){
 plot(year, dat[, i], type="l", ylab="", main=mains[i])
}
plot(year, residual, type="l", ylab="", main="residuals")
acf(residual, main="ACF of residuals", xlab="lag")

options(digits=6)
po.test(dat[, c(3, 1, 2, 4, 5)])

options(digits=3)
summary(ca.jo(dat))

library(fEcofin)
library(urca)
x = midcapD.ts[,2:11]
prices = exp(apply(x, 2, cumsum))
options(digits=3)
fit = summary(ca.jo(prices))
fit
loading = fit@W
alpha = fit@V %*% diag(fit@lambda) %*% solve(fit@V)
t(loading %*% alpha %*% t(x))


library("fEcofin")
library(urca)
# Short term rates
mk.maturity[2:11,]
fit = summary(ca.jo(mk.zero2[,2:11]))
fit

# Assuming:
#   A hedge fund owns a P0 portfolio
#   Ps amount is owned by the hedge fund, and P0 - Ps is borrowed.
#   At the end of the day, 
#     the fund will liquid the portfolio when its value fall bellow P0 - Ps
#     or when its value is above Pe
#   The portfolio was selected using cointegration analysis and its value is an AR(1) stationary process
#
# From above assumption:
#   (Pt - mu) = phi * (Pt_1 - mu) + St
#   Where:
#     Pt is the value of the portfolio at end of day t
#     mu is the mean of the portfolio
#     Pt_i is the value of the portfolio at end of day t - i
#     St is white noise with standard deviation sigma
#     phi < 1
#   Pt is a random normal variable with:
#    mean(Mut):                (P0 - mu) * phi^t + mu
#    standard deviation(Sdt):  sqrt((phi^(2*t) - 1) / (phi^2 - 1)) * sigma
#
# Let:
#    Lt: the possibility that the portfolio will be liquidated at the end of day t.
#    Nt: the possibility that the portfolio will be liquidated for a loss at the end of day t
#    Et: the possibility that the portfolio will be liquidated for a profit at the end of day t
#    Qet: the possibility that Pt >= Pe
#    Qlt: the possibility that Pt <= P0 - Ps
#    Rt: the condional expected return when the position is liquidated at the end of day t.
#    Slt = (1 - Lt)* (1 - Lt_1) * .... * (1 - L1)
# Then:
#    Nt = Qlt * Slt_1
#    Et = Qet * Slt_1
#    Lt = Nt + Et when t >= 1
#         0       when t == 0
#
# Quest1: What is the expected profit?
# Quest2: What is the probability that the hedge fund will need to liquidate for a loss?
# Quest3: What is the expected waiting time until the portfolio is liquidated?
# Quest4: What is the expected yearly return on the $50,000 investment?

# Approach 1: analitical way

P0 = 1e6
Ps = 5e4
Pe = 1.02e6
mu = 1.03e6
sigma = 5e3
phi = 0.99

# Let's do for 1000 days
days = 1000
Nt = rep(0, days)
Et = rep(0, days)
Lt = rep(0, days)
Rt = rep(0, days)
Log_Slt = rep(0, days)
Mu1 =  (P0 - mu) * phi + mu
Nt[1] = pnorm(P0 - Ps, mean=Mu1, sd=sigma)
Et[1] = 1 - pnorm(Pe, mean=Mu1, sd=sigma)
Lt[1] = Nt[1] + Et[1]
Log_Slt[1] = log(1 - Lt[1])
my_func = function(mu, x, sd){
  exp((-((x - mu) / sd)^2) / 2)
}
Rt[1] = Mu1 + sigma*(my_func(Mu1, Pe, sigma) - my_func(Mut, P0 - Ps, sigma)) / (sqrt(2 * pi) * (Nt[1] + Et[1]))

for(t in 2:days){
  Mut = (phi^t) * (P0 - mu) + mu
  Sdt = sqrt((phi^(2*t) - 1) / (phi^2 - 1)) * sigma
  Qlt = pnorm(P0 - Ps, mean=Mut, sd=Sdt)
  Qet = 1 - pnorm(Pe, mean=Mut, sd=Sdt)
  exp_log_Slt = exp(Log_Slt[t - 1])
  Nt[t] = Qlt * exp_log_Slt
  Et[t] = Qet * exp_log_Slt
  Lt[t] = Nt[t] + Et[t]
  Log_Slt[t] = Log_Slt[t - 1] + log(1 - Qlt - Qet)
  Rt[t] = Mut + Sdt*(my_func(Mut, Pe, Sdt) - my_func(Mut, P0 - Ps, Sdt)) / (sqrt(2 * pi) * (Qlt + Qet)) - P0
}

#Quest1:
sum(Lt*Rt)
# [1] 27859.56

#Quest2:
cumsum(Nt)[days]
# [1] 0.004557307

#Quest3:
sum((1:days)*Lt)
# [1] 11.78258

#Quest4:
sum(Lt*Rt / (Ps * (1:days))) * 253
# [1] 14.47078

# Approach 2: simulation
num_iterations = 100000
profits = rep(0, num_iterations)
expected_days = rep(0, num_iterations)
prices4 = c()
for(i in 1:num_iterations){
  price = P0
  t = 0
  while(price < Pe & price > P0 - Ps){
    t = t + 1
    Mut = (phi^t) * (P0 - mu) + mu
    Sdt = sqrt((phi^(2*t) - 1) / (phi^2 - 1)) * sigma
    price = rnorm(1, mean=Mut, sd=Sdt)
#   tried the following:
#     price = phi * (price - mu) + rnorm(1, mean=0, sd=sigma) + mu
#   The above code is biased against the analitical way.
#   The rationale here is(I guess):
#     price calculated in above manner depends on the fact that price in the previous iteration
#     is in the range [P0-Ps, Pe], this depenency makes it different from the analitical way.
# This is to check if the prices is normally distributed when the biased way
# is used.
    if(t == 4){
      prices4 = c(prices4, price - P0)
    }
  }
  profits[i] = price - P0
  expected_days[i] = t
}
# Why it's not normally distributed?
qqnorm(prices4)

#Quest1:
mean(profits)
# [1] 27792.19

#Quest2:
sum(profits < 0) / num_iterations
# [1] 0.00471

#Quest3:
mean(expected_days)
# [1] 11.77544

#Quest4:
mean(profits / (expected_days * Ps)) * 253
# [1] 14.26593

# price calculated in the following way is normally distributed with
#    mean = (phi^4) * (P0 - mu) + mu, sd = sqrt((phi^8 - 1)/(phi^2 - 1)) * sigma
num_iterations = 10000
prices = rep(0, num_iterations)
for(i in 1:num_iterations){
  P1 = phi * (P0 - mu) + rnorm(1, mean=0, sd=sigma) + mu
  P2 = phi * (P1 - mu) + rnorm(1, mean=0, sd=sigma) + mu
  P3 = phi * (P2 - mu) + rnorm(1, mean=0, sd=sigma) + mu
  P4 = phi * (P3 - mu) + rnorm(1, mean=0, sd=sigma) + mu
  prices[i] = P4 - P0
}

mean(prices)