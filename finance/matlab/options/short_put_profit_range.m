function [asset_prices, profits]=short_put_profit_range(strike_price, option_price, asset_price_start, asset_price_end, count)
	asset_prices = asset_price_start : (asset_price_end - asset_price_start) / count : asset_price_end;
	[~, d2] = size(asset_prices);
	profits = zeros(1, count);
	for i=1:d2
		if asset_prices(i) > strike_price
			profits(i) = option_price;
		else
			profits(i) = asset_prices(i) - strike_price + option_price;
		end
	end
end