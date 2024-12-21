%% LECTURE 8 - Monte Carlo Simulation - Noé Debrois - 16/10/2024
% MC with variance reduction technique : Antithetic Variables (AV).
% Pricing EU Call Option under B&S by MC with AV.
clear
clc
%% Parameters :
% MC parameter :
Nsim = 1e6;
% Contract parameters :
T = 1; K = 100;
% Market parameters :
S0 = 96; r = 0.01; 
% Market CALIBRATED parameter :
sigma = 0.6; % This parameter has to be calibrated.
%
%% Simulation :
% Source of randomness :
x1 = randn(Nsim, 1); % x1 ~ N(0, 1).
x2 = -x1; % Antithetic variable : perfectly negative-correlated RVs w/ x1.
% IMPLEMENTATION with U RV and icdf:
%U1 = rand(Nsim, 1);
%U2 = 1 - U1;
%x1 = icdf('normal', U1, 0, 1);
%x2 = icdf('normal', U2, 0, 1);

% B&S framework :
ST1 = S0 * exp((r - sigma^2 / 2) * T + sigma * sqrt(T) * x1);
ST2 = S0 * exp((r - sigma^2 / 2) * T + sigma * sqrt(T) * x2);

% Computation of the discounted payoffs :
discpayoff1 = exp(-r * T) * max(ST1 - K, 0);
discpayoff2 = exp(-r * T) * max(ST2 - K, 0);

% Correlations (FYI) :
corr(x1, x2) % This should be perfectly negative-correlated.
%corr(ST1, ST2)
%corr(discpayoff1, discpayoff2)
%
%% Pricing :
% Compute the MC price :
MC_Price = mean(1/2 * (discpayoff1 + discpayoff2));
MC_Price
% Get the exact price (B&S) to compare :
EXACT_Price = blsprice(S0, K, r, T, sigma);
EXACT_Price

% What if we don't have an exact price but we want to know the quality of
% the MC ? Let's use a Confidence Interval (CI) :
discpayoff = 1/2 * (discpayoff1 + discpayoff2);
[MC_Price, ~, CI] = normfit(discpayoff);
CI % By fixing a certain level of confidence, we can be sure on the i^th 
% decimal of the estimation.
%
%% PS :
% We just need to change the simulation of ST  if we want to change the
% model. The rest stays the same and the MC still works.
