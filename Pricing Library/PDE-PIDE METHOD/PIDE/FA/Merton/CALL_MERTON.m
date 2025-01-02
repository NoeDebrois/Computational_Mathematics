% ----------------------------------------------------------------------- %
clear; close all;
% Price an European Call option under the Merton model
% PIDE -> logprice, IMPLICIT scheme (operator splitting)
% ----------------------------------------------------------------------- %

%% Parameters & Grids (space + time):
% Option parameters:
S0 = 1; K = 1; T = 0.5; r = 0.1 / 100; sigma = 0.4;

% Merton model parameters:
lambda = 2; muJ = -0.01; deltaJ = 0.2;

% Time grid:
M = 250; dt = T / M;

% Space grid:
N = 1000;
Smin = 0.3 * S0; Smax = 3 * S0; 
xmin = log(Smin / S0); xmax = log(Smax / S0);
dx = (xmax - xmin) / N; 
x = xmin + (0:N)' * dx; % Notice "'" !

% CHANGE IF PUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the Lévy integral term, in the interpolation helper function:
UP_BC = @(t,x) S0 * exp(x) - K*exp(-r * t); % v(x,t) for x>xmax
LO_BC = @(t,x) 0;                           % v(x,t) for x<xmin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Integral truncation, quadrature points:
% Lévy Measure "nu" of Merton model:
nu = @(y)lambda * exp(-(y - muJ) .^2 / (2 * deltaJ^2)) / (deltaJ * sqrt(2 * pi));

% ----------------------------------------------------------------------- %
% Integral domain truncation:
tol = 1e-8;

% Set lb as "the y such that nu(y) < tol on the left": 
tmin = -0.5;
while (abs(nu(tmin)) > tol)
    tmin = tmin - 0.5;
end

% Set ub as "the y such that nu(y) < tol on the right": 
tmax = 0.5;
while (abs(nu(tmax)) > tol)
    tmax = tmax + 0.5;
end

% Plot of the kernel, i.e, the Lévy Measure between lb and ub:
figure
tt = linspace(tmin, tmax, 2 * N);
plot(tt, nu(tt));
title("Lévy Measure between l_b & u_b (i.e the 'kernel')", FontSize=15);
% ----------------------------------------------------------------------- %

%% Implicit matrices (A, B, C):
Mat = spalloc(N + 1, N + 1, 3 * (N - 1) + 2);
Mat(1, 1) = 1; Mat(end, end) = 1;

A = -(r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2);
B = -1 / dt - sigma^2 / (dx^2) - r;
C = (r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2);
for i=2:N
    Mat(i, [i-1 i i+1]) = [A B C];
end

%% Backward in time:
c = max(0, S0 * exp(x) - K); % Payoff @ maturity [CHANGE IF PUT]
rhs = zeros(N + 1, 1);
for j=M:-1:1 
    % CHANGE IF PUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rhs(1) = 0; % BC @ xmin
    rhs(end) = S0 * exp(xmax) - K * exp(-r * (T - (j - 1) * dt)); % BC @ xmax
    % We could also use UP_BC/LO_BC fcts with T - (j - 1) * dt & xmax/xmin
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Lévy integral term computation:
    ub_BC = @(x) UP_BC(T - j * dt, x);
    lb_BC = @(x) LO_BC(T - j * dt, x);
    I = INTEGRAL_LEVY(x, c, tt, nu, lb_BC, ub_BC);

    rhs(2:end-1) = -1 / dt * c(2:end-1) - I;
    c = Mat\rhs;
end

%% Results:
% Plot of the price of the option vs a grid of spot prices:
figure
plot(S0 * exp(x), c);
title("EU (call or put) option price VS a grid of spot prices", FontSize=15);

% Interpolate at the specific spot price we entered:
price = interp1(x, c, 0, 'spline')
% Check with Carr-Madan price:
price_cm = FFT_CM_Call_Merton(K, [sigma lambda muJ deltaJ], T, r, S0)