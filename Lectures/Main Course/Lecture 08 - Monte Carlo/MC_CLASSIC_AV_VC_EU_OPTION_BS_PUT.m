function Price = MC_CLASSIC_AV_VC_EU_OPTION_BS_PUT(S0, K, r, T, sigma, Nsim, flag, Nsim2)
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

nuT = (r - sigma^2 / 2) * T;
sigmaT = sigma * sqrt(T); % sigma * sqrt(T) * N(0,1) =(d) N(0,sigma^2 * T)
temp = randn(Nsim, 1);

if flag == 0 % NO VARIANCE REDUCTION
    X = S0 * exp(nuT + sigmaT * temp);
    DiscPayoffs = exp(-r * T) * max(K - X, 0); % Put
    Price = sum(DiscPayoffs) / Nsim;

elseif flag == 1 % ANTITHETIC VARIABLES
    X1 = S0 * exp(nuT + sigmaT * temp);
    X2 = S0 * exp(nuT + sigmaT * (-temp));
    Payoff1 = max(K - X1, 0);
    Payoff2 = max(K - X2, 0);
    DiscPayoffs = exp(-r * T) * 1/2 * (Payoff1 + Payoff2);
    Price = DiscPayoffs / Nsim;

elseif flag == 2 % CONTROL VARIABLES
    % 1 - Compute parameters (alpha & E[f])
    temp2 = randn(Nsim2, 1);
    % -------- f & g -------- %
    f = S0 * exp(nuT + sigmaT * temp2); % Choice of f: underlying asset
    g = exp(-r * T) * max(K - f, 0);   % g for the estimation of alpha
    % -------- alpha -------- %
    Ef = S0 * exp(r * T); % Expectation of f
    Varf = S0^2 * exp(2 * r * T) * (exp(T * sigma^2) - 1); % Variance of f
    Covfg = cov(f, g);
    alpha = - Covfg(1, 2) / Varf; % cf formula in the notes

    % 2 - Compute price
    % -------- f & g -------- %
    New_f = S0 * exp(nuT + sigmaT * temp); % f but for the "true" MC
    New_g = exp(-r * T) * max(K - New_f, 0); % Discounted payoffs
    DiscPayoffs = New_g + alpha * (New_f - Ef);
    Price = sum(DiscPayoffs) / Nsim;
end

% To get a confidence interval:
[Price, ~, CI] = normfit(DiscPayoffs);
disp("Confidence interval:")
disp(CI);

% For comparison with the closed formula:
[~, BSPrice] = blsprice(S0, K, r, T, sigma);
disp("BS price:")
disp(BSPrice)
disp("Difference:")
disp(Price - BSPrice);