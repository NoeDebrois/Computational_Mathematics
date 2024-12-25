function [call_price, call_CI_price, put_price, put_CI_price] = MC_MERTON(S0, K, r, T, sigma, lambda, muJ, deltaJ, M, Nsim, flag, Nsim2)
% S0    : price of the underlying asset
% K     : strike
% r     : risk-free interest rate 
% T     : maturity 
% sigma : volatility
% Nsim  : number of Monte Carlo simulations
% flag  : type of variance reduction considered
%  --> flag = 0 : no variance reduction
%  --> flag = 1 : antithetic variables 
%  --> flag = 2 : control variables
% Nsim2 : number of MC simulations to compute the alpha (only if flag = 2)

if flag == 0 % NO VARIANCE REDUCTION
    [~, ST] = MERTON_ASSET(S0, r, T, sigma, lambda, muJ, deltaJ, M, Nsim);
    DiscPayoffs = exp(-r * T) * max(ST - K, 0); % Call

elseif flag == 1 % ANTITHETIC VARIABLES
    [~, ~, ST, STAV] = MERTON_ASSET_AV(S0, r, T, sigma, lambda, muJ, deltaJ, M, Nsim);
    Payoff1 = max(ST - K, 0);
    Payoff2 = max(STAV - K, 0);
    DiscPayoffs = exp(-r * T) * 1/2 * (Payoff1 + Payoff2);

elseif flag == 2 % CONTROL VARIABLES
    % 1 - Compute parameters (alpha & E[f])
    % -------- f & g -------- %
    f = MERTON_ASSET(S0, r, T, sigma, lambda, muJ, deltaJ, M, Nsim2);
    g = exp(-r * T) * max(f - K, 0); % f & g for the estimation of alpha
    % -------- alpha -------- %
    Ef = S0 * exp(r * T); % f = ST, so under Q: E[f(Z)] = S0 * exp(rT)
    Covfg = cov(f, g);
    alpha = - Covfg(1, 2) / Covfg(1, 1); % cf formula in the notes

    % 2 - Compute price
    % -------- f & g -------- %
    [~, New_f] = MERTON_ASSET(S0, r, T, sigma, lambda, muJ, deltaJ, M, Nsim);
    New_g = exp(-r * T) * max(New_f - K, 0); % f & g for the price
    DiscPayoffs = New_g + alpha * (New_f - Ef);
end

% To get a confidence interval:
[call_price, ~, call_CI_price] = normfit(DiscPayoffs)

% Compute put price and asymptotic CI via Call-Put Parity:
put_parity = @(call_p) call_p - S0 + K*exp(-r*T);
put_price = put_parity(call_price);
put_CI_price = put_parity(call_CI_price);

% For comparison with Carr-Madan price:
params.sigma = sigma;
params.mu = muJ;
params.delta = deltaJ;
params.lambdaK = lambda;
CHECK_CM_CALL_PRICE = MERTON_CARR_MADAN(K, params, T, r, S0);
disp("[CHECK] Carr-Madan call price:");
disp(CHECK_CM_CALL_PRICE);