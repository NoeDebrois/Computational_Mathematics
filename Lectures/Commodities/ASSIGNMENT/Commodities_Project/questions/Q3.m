%% Q3 - Model Calibration on 2026 prices

clc
format default

addpath("data\")
addpath("utility\")
addpath("plot\")

% Initial time
t0 = data.t0;

% Initial forward price
F0 = data.fwd(data.fwd.Name == "2026", :).Price;

% Volatility surface
T = t0 + calyears(ceil(data.volsurface_2026.tenors));
T = floor(data.volsurface_2026.tenors.*ndays(t0, T));
T = t0 + T;
K = data.volsurface_2026.strikes';
volatility = data.volsurface_2026.vols;

% Discount rates 
[~, r] = discounts(data.discounts(1,:), data.discounts(2,:), t0, T);
[~, ttm] = ndays(t0,T);

% Market prices
price_BS = BS(volatility, F0, K, r, ttm);

%% Calibrate the HJM on the mrkt prices 
% FFT paramters
A = 1000; n = 16;
% Calibration parameters
model = 'NIG';
time_dep = 'none'; 

% Initial parameters
x0 = [ 0.5,   0.1,     0.5,    0.1];
lb = [0.01,    -2,    0.01,     -2];        
ub = [   1,    +2,       5,     +2];

% Allow upsilon to be time dependent
switch(time_dep)
    case('none')
    case('upsilon')
        x0 = [x0, repmat(x0(end), 1, 6)];
        lb = [lb, repmat(lb(end), 1, 6)];
        ub = [ub, repmat(ub(end), 1, 6)]; 
    otherwise
        return
end

% Calibration
options =  optimoptions("lsqnonlin", "Display", "off", "MaxFunctionEvaluations", 1e5, "MaxIterations", 1e3);
[params, ~] = lsqnonlin( @(params) ...
                  price_err(volatility, F0, K, r, ttm, A, n, params, model, time_dep), ...
                  x0, lb, ub, options);
error = mean(abs(price_err(volatility, F0, K, r, ttm, A, n, params, model, time_dep)), 'all');

%% Plot
% Plot price surfaces
price_HJM = HJM(A, n, F0, r, ttm, K, params, model, time_dep);
plot_price_3d(K, ttm, price_BS, price_HJM);

% Plot implied volatility surfaces
volatility_HJM = implied_vol(price_HJM, F0, K, r, ttm);
plot_vol_3d(K, ttm, volatility, volatility_HJM);

%% Save calibrated parameters
data.parameters.Q3 = params;

%% Print results
fprintf("Q3 - Calibrated parameters: \n")
fprintf("Sigma = %.4f \n", params(1))
fprintf("Theta = %.4f \n", params(2))
fprintf("Kappa = %.4f \n", params(3))
fprintf("Upsilon = %.4f \n", params(4))
fprintf("\nMAE = %.2f \n", error)
