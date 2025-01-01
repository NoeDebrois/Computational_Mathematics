%% Q6 - Option pricing 2026 & 2028

clc
format default
rng(0)

addpath('data\')
addpath('utility\')
addpath('plot\')

% load calibrated parameters
params = data.parameters.Q5;

% Initial time
t0 = data.t0;

% Initial forward prices
F0_26 = data.fwd(data.fwd.Name == "2026", :).Price;
F0_28= data.fwd(data.fwd.Name == "2028", :).Price;
K = 300;

% Maturity
T = t0 + calyears(1);

% Discount rates 
[~, r] = discounts(data.discounts(1,:), data.discounts(2,:), t0, T);
[~, ttm] = ndays(t0,T);

% Option price
N = 1e5;
M = floor(12*ttm);

% MC
method = 'MC';
[price_MC, CI_MC] = option_26_28(N, M, F0_26, F0_28, K, r, ttm, params, method);

% AV
method = 'AV';
[price_AV, CI_AV] = option_26_28(N, M, F0_26, F0_28, K, r, ttm, params, method);

%% Print results
fprintf("\nQ6 - Option price: \n")
fprintf('Price MC: %.2f CI: [%.2f, %.2f] \n', price_MC, CI_MC);
fprintf('Price AV: %.2f CI: [%.2f, %.2f] \n', price_AV, CI_AV);
