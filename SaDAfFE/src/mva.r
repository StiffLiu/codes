mva = function(ar, samples){
  total = length(ar)
  len = min(total, samples)
  beg = cumsum(ar[1:len]) / (1:len)
  if (total > len){
	count = total - samples
	end = rep(0, count)
	for (i in 1:count) end[i] = mean(ar[i:(i + samples)])
	beg = c(beg, end)
  }
  beg
}
dat = read.csv("F:\\315.csv")
prices = dat[[1]]
plot(prices, type='l')
lines(mva(prices, 1000), type='l', col='red')
lines(mva(prices, 2000), type='l', col='green')
lines(mva(prices, 5000), type='l', col='blue')
lines(mva(prices, 10000), type='l', col='yellow')