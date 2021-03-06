addpath('options');
strike_price = 100;
option_price = 5;
asset_start = strike_price - 2 * option_price;
asset_end = strike_price + 2 * option_price;
count = 100;
[sc_r,short_call] = short_call_profit_range(strike_price, option_price, asset_start, asset_end, count);
[lc_r,long_call] = long_call_profit_range(strike_price, option_price, asset_start, asset_end, count);
[sp_r,short_put] = short_put_profit_range(strike_price, option_price, asset_start, asset_end, count);
[lp_r,long_put] = long_put_profit_range(strike_price, option_price, asset_start, asset_end, count);
plot(sc_r, short_call, lc_r, long_call, sp_r, short_put, lp_r, long_put);
xlim([asset_start - 1, asset_end + 1]);
ylim([-option_price * 1.1, option_price * 1.1]);
legend('Short Call', 'Long Call', 'Short Put', 'Long Put');

% subplot(1, 2, 1);
% plot(sc_r, short_call, lc_r, long_call)
% xlim([asset_start - 1, asset_end + 1]);
% ylim([-option_price * 1.1, option_price * 1.1]);
% legend('Short', 'Long');
% title('Call option');
% 
% subplot(1,2,2)
% plot(sp_r, short_put, lp_r, long_put)
% xlim([asset_start - 1, asset_end + 1]);
% ylim([-option_price * 1.1, option_price * 1.1]);
% legend('Short', 'Long');
% title('Put option');