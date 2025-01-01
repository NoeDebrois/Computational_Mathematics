%% Q5 - Model Calibration on 2026 & 2028 prices

clc
format default

addpath("data\")
addpath("utility\")
addpath("plot\")

% Initial time
t0 = data.t0;

% Initial forward prices
F0_26 = data.fwd(data.fwd.Name == "2026", :).Price;
F0_28= data.fwd(data.fwd.Name == "2028", :).Price;

% Volatility surface 2026
T_26 = t0 + calyears(ceil(data.volsurface_2026.tenors));
T_26 = floor(data.volsurface_2026.tenors.*ndays(t0, T_26));
T_26 = t0 + T_26;
K_26 = data.volsurface_2026.strikes';
volatility_26 = data.volsurface_2026.vols;

% Volatility surface 2026
T_28 = t0 + calyears(ceil(data.volsurface_2028.tenors));
T_28 = floor(data.volsurface_2028.tenors.*ndays(t0, T_28(1)));
T_28 = t0 + T_28;
K_28 = data.volsurface_2028.strikes';
volatility_28 = data.volsurface_2028.vols;

% Discount rates 2026
[~, r_26] = discounts(data.discounts(1,:), data.discounts(2,:), t0, T_26);
[~, ttm_26] = ndays(t0,T_26);

% Discount rates 2026
[~, r_28] = discounts(data.discounts(1,:), data.discounts(2,:), t0, T_28);
[~, ttm_28] = ndays(t0,T_28);

% Market prices 2026
price_BS_26 = BS(volatility_26, F0_26, K_26, r_26, ttm_26);

% Market prices 2028
price_BS_28 = BS(volatility_28, F0_28, K_28, r_28, ttm_28);

%% Calibrate the HJM on the mrkt prices 
% FFT paramters
A = 1000; n = 16;
% Calibration parameters
model = 'NIG';
time_dep = 'none'; 

x0 = [ 0.5,   0.1,     0.5,    0.1];
lb = [0.01,    -2,    0.01,     -2];        
ub = [   1,    +2,       5,     +2];

options =  optimoptions("lsqnonlin", "Display", "off", "MaxFunctionEvaluations", 1e5, "MaxIterations", 1e3);
[params, ~] = lsqnonlin( @(params) ...
                  [price_err(volatility_26, F0_26, K_26, r_26, ttm_26, A, n, params, model, time_dep); ...
                   price_err(volatility_28, F0_28, K_28, r_28, ttm_28, A, n, params, model, time_dep)], ...
                  x0, lb, ub, options);
error = mean(abs([price_err(volatility_26, F0_26, K_26, r_26, ttm_26, A, n, params, model, time_dep); ...
                  price_err(volatility_28, F0_28, K_28, r_28, ttm_28, A, n, params, model, time_dep)]),  'all');

%% Plot
% Plot price surfaces 2026
price_HJM_26 = HJM(A, n, F0_26, r_26, ttm_26, K_26, params, model, time_dep);
plot_price_3d(K_26, ttm_26, price_BS_26, price_HJM_26);

% Plot implied volatility surfaces 2026
volatility_HJM_26 = implied_vol(price_HJM_26, F0_26, K_26, r_26, ttm_26);
plot_vol_3d(K_26, ttm_26, volatility_26, volatility_HJM_26);

% Plot price surfaces 2026
price_HJM_28 = HJM(A, n, F0_28, r_28, ttm_28, K_28, params, model, time_dep);
plot_price_3d(K_28, ttm_28, price_BS_28, price_HJM_28);

% Plot implied volatility surfaces 2026
volatility_HJM_28 = implied_vol(price_HJM_28, F0_28, K_28, r_28, ttm_28);
plot_vol_3d(K_28, ttm_28, volatility_28, volatility_HJM_28);

%% Save parameters in data
data.parameters.Q5 = params;

%% Print results
fprintf("Q5 - Calibrated parameters: \n")
fprintf("Sigma = %.4f \n", params(1))
fprintf("Theta = %.4f \n", params(2))
fprintf("Kappa = %.4f \n", params(3))
fprintf("Upsilon = %.4f \n", params(4))
fprintf("\nMAE = %.2f \n", error)
