num_iteration = 1e5
days = 100
trading_days_per_year = 253
current_price = 1e6
log_loss = log(0.95)
log_gain = log(1.1)
mean_log_ret_per_day = 0.05 / trading_days_per_year
var_log_ret_per_day = 0.23 / sqrt(trading_days_per_year)
margin_price = 9.5e5

gain_distribution = rep(0, num_iteration)
loss_distribution = rep(0, num_iteration)
profits = rep(NA, num_iteration)
returns = rep(NA, num_iteration)

for (i in 1:num_iteration){
  #generate accumulated random log returns for the following _days_ days
  a_l_r = cumsum(rnorm(days, mean=mean_log_ret_per_day, sd=var_log_ret_per_day))
  #get the profit when the stock will be sold
  log_ret = c(a_l_r[a_l_r <= log_loss | a_l_r >= log_gain], a_l_r[days])[1]
  gain_distribution[i] = as.numeric(log_ret >= log_gain)
  loss_distribution[i] = as.numeric(log_ret <= 0)
  profits[i] = exp(log_ret) * current_price
  if (profits[i] > margin_price) returns[i] = log((profits[i] - margin_price) / (current_price - margin_price)) / i
}

#profits = profits - current_price
returns = returns[!is.na(returns)]
mean(gain_distribution)
mean(loss_distribution)
mean(profits)
mean(returns)