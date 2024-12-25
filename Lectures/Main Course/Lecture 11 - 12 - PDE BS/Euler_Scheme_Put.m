%% LECTURE 11 - Explicit Euler for EU Put - Noé Debrois - 29/10/2024
% Script to price a EU Put Option, under B&S model, using log-price PDE.
% We implement Finite Difference Method using EXPLICIT EULER SCHEME.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About EXPLICIT EULER SCHEME :
% We make an error at each iteration. The Explicit Euler Scheme is instable
% because the errors cumulate each other. And we do not know exactly how
% large M shoud be. It is "M >> N - conditionally stable" : try and error.
% Explicit Euler is VERY SIMPLE but CONDITIONALLY STABLE so let's take a
% look at "Implicit Euler".
clear
close all
%
%% 1. Parameters & Grids :
% Model parameters :
S0 = 1; K = 0.95; r = 0.001; sigma = 0.5; T = 1;

% Space grid :
N = 100; % Number of points in the space dimension.
Smin = 0.1 * S0; % Found by experience : we assume that the spot price has
Smax = 3 * S0;   % a negligible probability of falling below 10% of S0, and
                 % of rising above 3 * S0.
xmin = log(Smin/S0); % x lives in [xmin, xmax]=[log(Smin/S0),log(Smax/S0)].
xmax = log(Smax/S0);
dx = (xmax - xmin) / N; % Space grid step.
x = xmin + (0:N) * dx;  % From x0 = xmin to xN = xmax (in total : N+1 pts).
% x = linspace(xmin, xmax, N+1);

% Time grid :
M = 1000;   % Number of points in the time dimension.
dt = T / M; % Time grid step.
t = (0:M) * dt; % t lives in [0, T], divided in M subintervals : M+1 pts.
%
%% 2. Explicit Euler - Backward in time loop :
% Vector of price (for all x) that we will update backwardly from T to 0 
% (i.e from tM to t0) :
c = max(K - S0 * exp(x'), 0);

cnew = zeros(size(c)); % New price vector that we will compute using prices
% from c (in cnew : prices at time t{j-1} ; in c : prices at time t_{j}).

for j=M:-1:1 % BACKWARD from t_M to t_0.
    % We know c t_j -> we compute cnew t_{j-1}.
    %
    % Boundary condition for PUT : "v(t, xmin) = K * e^[-r(T-t)] - Smin, 
    % for all t" :
    cnew(1) = K * exp(-r * (T - (j - 1) * dt)) - Smin;
    %
    for i=2:N
        % Explicit Euler Scheme : by moving V_{j-1, i} (UNKNOWN) on the left :
        % V_{j-1, i} = -A * dt * V_{j, i-1} - B * dt * V_{j, i} - C * dt *
        % V_{j, i+1}, 
        %
        % Where V_{j, i-1}, V_{j, i} & V_{j, i+1} are known.
        %
        % And :
        %       A = (r - sigma^2 / 2) / (2 * dx) - (sigma^2 / 2) / dx^2 ;
        %       B = - 1 / dt + sigma^2 / dx^2 + r ;
        %       C = - (r - sigma^2 / 2) / (2 * dx) - (sigma^2 / 2) / dx^2.
        %
        cnew(i) = c(i) + dt * (...
            (r - sigma^2 / 2) * (c(i + 1) - c(i - 1)) / (2 * dx) +...
            sigma^2 / 2 * (c(i + 1) - 2 * c(i) + c(i - 1)) / dx^2 +...
            -r * c(i));
    end
    %
    % Boundary condition for PUT : "v(t, xmax) = 0, for all t" :
    cnew(N + 1) = 0;

    % Update the price @ t_j to the new price @ t_{j-1} (BACKWARD) :
    c = cnew;
end
%
%% Plot "put price at t=0" VS "S at t=0" (for all S in [Smin, Smax]) :
% S at t = 0 represents the spot price (so, now). S0 (cf "Model
% parameters") is one example of possible spot price.
figure
plot(S0 * exp(x), c); 
title('Put price at t=0 vs (all possible) Spot Price S_0');
xlabel('S at t = 0 (spot price)');
ylabel('Put price at t = 0');
%
%% Compute the price at t = 0 :
price = interp1(x, c, 0, 'spline') % Interpolation of c at x = 0, i.e S = S0
% (S0 being the one we the one we initialized in "Model parameters").
[~, price_ex] = blsprice(S0, K, r, T, sigma) % For comparison : B&S price.