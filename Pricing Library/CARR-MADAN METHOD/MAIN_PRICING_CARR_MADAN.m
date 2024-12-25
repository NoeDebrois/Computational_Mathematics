clear; close all;
addpath("CHAR_EXP/");

%% Plain Vanilla EU Call Option parameters:
Strike = [100];
S0 = 100;    % Initial Spot Price
param.T = 2; % Maturity

param.rf = 0.02; % Risk-free rate
param.q = 0.00;  % Dividend yield

%% Model parameters:
distr = 5; 

%% Pricing:

if distr == 1 % B&S
    % sigma = volatility of the model
    param.sigma = 0.2;
    BS_CARR_MADAN(Strike, param, param.T, param.rf, S0)

elseif distr == 2 % Merton
    % sigma   = BM volatility (as in GBM)
    % mu      = jump drift (meansize of jumps)
    % delta   = jump vol (std of jumps)
    % lambdaK = jump time intensity
    param.sigma   = 0.2;
    param.mu      = 0.01;
    param.delta   = 0.2;
    param.lambdaK = 2;
    MERTON_CARR_MADAN(Strike, param, param.T, param.rf, S0)

elseif distr == 3 % Kou
    % sigma   = BM volatility (as in GBM)
    % p       = prob of positive jump
    % lambdap = pos jump intensity
    % lambdam = neg jump intensity
    % lambdaK = jump time intensity
    param.sigma   = 0.2;
    param.p       = 0.6;
    param.lambdap = 10;
    param.lambdam = 3;
    param.lambdaK = 2;
    KOU_CARR_MADAN(Strike, param, param.T, param.rf, S0)

elseif distr == 4 % NIG
    % sigma = volatility of the BM
    % theta = drift of the BM
    % kNIG  = variance of the subordinator
    param.sigma = 0.2;
    param.theta = 0.05;
    param.kNIG  = 1; % kNIG->0^+ : like B&S.
    NIG_CARR_MADAN(Strike, param, param.T, param.rf, S0)

elseif distr == 5 % VG
    % sigma = volatility of the BM
    % theta = drift of the BM
    % k     = variance of the subordinator
    param.sigma = 0.6;
    param.theta = 0.2;
    param.kVG   = 0.5;
    VG_CARR_MADAN(Strike, param, param.T, param.rf, S0)

%elseif param.distr == 6
%elseif param.distr == 7
end