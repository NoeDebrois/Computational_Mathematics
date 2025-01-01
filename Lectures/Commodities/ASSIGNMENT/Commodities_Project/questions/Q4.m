%% Q4 - Option pricing 2026

clc
format default
rng(0)

addpath('data\')
addpath('utility\')
addpath('plot\')

% load calibrated parameters from Q3
params = data.parameters.Q3;

% Initial time
t0 = data.t0;

% Initial forward price
F0 = data.fwd(data.fwd.Name == "2026", :).Price;
K = 300;

% Maturity
T = t0 + calyears(1);

% Discount rates 
[~, r] = discounts(data.discounts(1,:), data.discounts(2,:), t0, T);
[~, ttm] = ndays(t0,T);

% Nig simulation parameters
N = 1e5;
M = floor(12*ttm);

% MC pricing
method = 'MC';
[price_MC, CI_MC] = option_26(N, M, F0, K, r, ttm, params, method);

% AV pricing (AV only applied to IG simulation)
method = 'AV';
[price_AV, CI_AV] = option_26(N, M, F0, K, r, ttm, params, method);

%% Print results
fprintf("\nQ4 - Option price: \n")
fprintf('Price MC: %.2f CI: [%.2f, %.2f] \n', price_MC, CI_MC);
fprintf('Price AV: %.2f CI: [%.2f, %.2f] \n', price_AV, CI_AV);
