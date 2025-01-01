function [price, CI] = option_26_28(N, M, F0_26, F0_28, K, r, T, params, method)

%--------------------------------------------------------------------------

% MC pricing for the option presented in Q6 (optional)

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

X_26 = zeros(N,M+1);
X_28 = zeros(N,M+1);
g_26 = randn(N,M);
g_28 = randn(N,M);

% Simulate from the inverse gaussian
U_26 = rand(N,M);
U_28 = rand(N,M);
dt = T/M;
mu = dt;
lambda = dt^2/k;
dS_26 = icdf("InverseGaussian", U_26, mu, lambda);
dS_28 = icdf("InverseGaussian", U_28, mu, lambda);
% Alternative sampling
%dS = IG(N, M, mu, lambda);

% Risk neutral drift
[~, drift, ~] = NIG_char_exp(params);
drift = r + drift;

% Simulate X according to NIG
for i = 1:M
    X_26(:,i+1) = X_26(:,i) + drift*dt + upsilon*(theta*dS_26(:,i) + sigma*sqrt(dS_26(:,i)).*g_26(:,i));
    X_28(:,i+1) = X_28(:,i) + drift*dt + upsilon*(theta*dS_28(:,i) + sigma*sqrt(dS_28(:,i)).*g_28(:,i));
end
F_26 = F0_26*exp(X_26);
F_28 = F0_28*exp(X_28);
[check_26,~,CI_check_26] = normfit(X_26(:,end)-r*T);
[check_28,~,CI_check_28] = normfit(X_28(:,end)-r*T);

% Compute the price
payoff = exp(-r*T)*max(max(F_26(:,end), F_28(:,end)) - K, 0);

% Simulate the anthitetic variable
if strcmp(method, 'AV')
    X_AV_26 = zeros(N,M+1);
    X_AV_28 = zeros(N,M+1);
    U_AV_26 = 1 - U_26;
    U_AV_28 = 1 - U_28;
    dS_AV_26 = icdf("InverseGaussian", U_AV_26, mu, lambda);
    dS_AV_28 = icdf("InverseGaussian", U_AV_28, mu, lambda);
    for i = 1:M
        X_AV_26(:,i+1) = X_AV_26(:,i) + drift*dt + upsilon*(theta*dS_AV_26(:,i) + sigma*sqrt(dS_AV_26(:,i)).*g_26(:,i));
        X_AV_28(:,i+1) = X_AV_28(:,i) + drift*dt + upsilon*(theta*dS_AV_28(:,i) + sigma*sqrt(dS_AV_28(:,i)).*g_28(:,i));
    end
    F_AV_26 = F0_26*exp(X_AV_26);
    F_AV_28 = F0_28*exp(X_AV_28);
    [check_26,~,CI_check_26] = normfit([X_26(:,end); X_AV_26(:,end)]-r*T);
    [check_28,~,CI_check_28] = normfit([X_28(:,end); X_AV_28(:,end)]-r*T);
    payoff_AV = exp(-r*T)*max(max(F_AV_26(:,end), F_AV_28(:,end)) - K, 0);
    payoff = (payoff + payoff_AV)/2;
end

% Normality check
fprintf('Normality check 2026: E[X_T - r*T] =  %.2f CI: [%.2f, %.2f] \n', check_26, CI_check_26);
fprintf('Normality check 2028: E[X_T - r*T] =  %.2f CI: [%.2f, %.2f] \n', check_28, CI_check_28);

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


