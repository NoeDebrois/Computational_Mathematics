function [driftless_psi, rn_drift, psi] = NIG_char_exp(params)

%--------------------------------------------------------------------------

% Compute NIG CharExp parametrically with r.n drift

% INPUT:
% params >> sigma, theta, kappa, upslion (from HJM)

% OUTPUT:
% drfitless_psi >> CharExp w/o drift
% rn_drift
% psi >> CharExp with drift

%--------------------------------------------------------------------------

% Parameters
sigma = params(1);
theta = params(2); kappa = params(3);
upsilon = params(4);

% Driftless CharExp
driftless_psi = @(u) (1 - sqrt(1 + u.^2*upsilon^2*sigma^2*kappa - ...
                      2*1i*theta*u*upsilon*kappa))./kappa;

% R.N condition
rn_drift = -driftless_psi(-1i);

% CharExp with R.N. condition
psi = @(u) 1i*u*rn_drift + driftless_psi(u);

end


