function [price, CI] = option_26(N, M, F0, K, r, T, params, method)

%--------------------------------------------------------------------------

% MC pricing for the option presented in Q4

% INPUT
% N >> Number of simulations
% M >> Number of time steps
% F0
% K
% r
% T
% params >> sigma, theta, kappa, upsilon

% OUTPUT 
% price
% CI

%--------------------------------------------------------------------------

% NIG parameters
sigma = params(1);
theta = params(2); k = params(3);
upsilon = params(4);

X = zeros(N,M+1);
g = randn(N,M);

% Simulate from the inverse gaussian
U = rand(N,M);
dt = T/M;
mu = dt;
lambda = dt^2/k;
dS = icdf("InverseGaussian", U, mu, lambda);
% Alternative sampling
%dS = IG(N, M, mu, lambda);

% Risk neutral drift
[~, drift, ~] = NIG_char_exp(params);
drift = r + drift;

% Simulate X according to NIG
for i = 1:M
    X(:,i+1) = X(:,i) + drift*dt + upsilon*(theta*dS(:,i) + sigma*sqrt(dS(:,i)).*g(:,i));
end
F = F0*exp(X);
[check,~,CI_check] = normfit(X(:,end)-r*T);

% Compute the price
payoff = exp(-r*T)*max(max(F, [], 2) - K, 0);

% Simulate the anthitetic variable
if strcmp(method, 'AV')
    X_AV = zeros(N,M+1);
    U_AV = 1 - U;
    dS_AV = icdf("InverseGaussian", U_AV, mu, lambda);
    for i = 1:M
        X_AV(:,i+1) = X_AV(:,i) + drift*dt + upsilon*(theta*dS_AV(:,i) + sigma*sqrt(dS_AV(:,i)).*g(:,i));
    end
    F_AV = F0*exp(X_AV);
    [check,~,CI_check] = normfit([X(:,end); X_AV(:,end)]-r*T);
    payoff_AV = exp(-r*T)*max(max(F_AV, [], 2) - K, 0);
    payoff = (payoff + payoff_AV)/2;
end

% Normality check
fprintf('Normality check: E[X_T - r*T] =  %.2f CI: [%.2f, %.2f] \n', check, CI_check);

% Pricing
[price, ~, CI] = normfit(payoff);

end

%--------------------------------------------------------------------------

% function x = IG(N, M, mu, lambda)
% 
% %------------------------------------------------------------------------
% 
% % Simulate from an IG(mu,lambda)
% 
% %------------------------------------------------------------------------
% 
% g = randn(N,M);
% g = g.^2;
% x = mu + (mu^2*g)/(2*lambda) - (mu/(2*lambda))*sqrt(4*mu*lambda*g + mu^2*g.^2);
% U = randn(N,M);
% idx = find(U > mu./(mu + x));
% x(idx) = mu^2./x(idx);
% 
% end

