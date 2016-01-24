function [asset_prices, profits] = long_put_profit(strike_price, option_price)
	[asset_prices, profits] = long_put_profit_range(strike_price, option_price,  strike_price - option_price, strike_price + 2 * option_price, 100);
end