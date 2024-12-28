%% LECTURE 11 - Theta Method for Up&Out Call - Noé Debrois - 29/10/2024
% Script to price a Up&Out Call Option, under B&S model, using price PDE.
% WARNING : NOT log-price PDE. cf Lecture notes, page 23 PDE.pdf.
% We implement Finite Difference Method using THETA METHOD SCHEME.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
%% 1. Parameters & Grids :
% Model & Option parameters :
S0 = 1; K = 0.95; r = 0.001; sigma = 0.5; T = 1;
U = 1.2; % Barrier.

% Method parameter :
theta = 0.5; 
% 0 : Implicit Euler, 0.5 : Crank-Nicholson, 1 : Explicit Euler

% Space grid :
N = 2000; % Number of points in the space dimension.
%
Smin = 0.1 * S0; % Found by experience : we assume that the spot price has 
% a negligible probability of falling below 10% of S0.
Smax = U; % The upper bound has to be THE BARRIER this time !
%
dS = (Smax - Smin) / N; % Space grid (price grid) step.
S = Smin + (0:N) * dS; % Price grid.

% Time grid :
M = 500;       % Number of points in the time dimension.
dt = T / M;     % Time grid step.
t = (0:M) * dt; % t lives in [0, T], divided in M subintervals : M+1 pts.
%
%% 2. Matrix Construction :
% We write the theta method in matrix form, cf my notebook about PDE.

% WARNING : with the PRICE PDE (contrary to logprice PDE), A, B & C depend
% on i. Hence the following functions to build A, B, C for each S_i :
A = @(S) (1 - theta) * dt * (-r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));
B = @(S) -1 + dt * (1 - theta) * (-sigma^2 * S^2 / (dS^2) - r);
C = @(S) dt * (1 - theta) * (r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));
% cf page 26 of notes PDE.pdf.

Mat = spalloc(N + 1, N + 1, 3 * (N - 1) + 2);
Mat(1, 1) = 1;
for i=2:N
    Mat(i, [i - 1, i, i + 1]) = [A(S(i)) B(S(i)) C(S(i))];
end
Mat(end, end) = 1;

% [Atilde, Btilde, Ctilde] (once again, with PRICE PDE, they depend on i) :
Ah = @(S) -(theta) * dt * (-r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));
Bh = @(S) -1 - dt * (theta) * (-sigma^2 * S^2 / (dS^2) - r);
Ch = @(S) -dt * (theta) * (r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));

Mat_rhs = spalloc(N + 1, N + 1, 3 * (N - 1));
for i=2:N
    Mat_rhs(i, [i - 1, i, i + 1])=[Ah(S(i)) Bh(S(i)) Ch(S(i))];
end
%
%% 3. Backward in time loop :
c = max(S' - K, 0) .* (S' < U) ; % @ Maturity : "v(T, x) = (ST - K)^+ * 1{S < U}". 

rhs = zeros(N + 1, 1); % BCj in my notes.

for j=M:-1:1
    % know c t_j -> compute cnew t_{j-1}
    rhs = Mat_rhs * c;
    rhs(1) = 0;   % BC at tj, for S=Smin -> CALL OPTION.
    rhs(end) = 0; % BC at tj, for S=Smax=U -> Hit the BARRIER.
    c = Mat \ rhs;
end
%
%% Plot "call price at t=0" VS "S at t=0" (for all S in [Smin, Smax]) :
figure
plot(S,c); 
title('U&O Call price at t=0 vs (all possible) Spot Price S_0'); 
xlabel('S at t = 0 (spot price)');
ylabel('U&O Call price at t = 0')
%
%% Compute the price at t = 0 :
price = interp1(S, c, S0, 'spline')

% RMK about apparent instability for Crank-Nicholson scheme (cf p30) :
% Payoff is not continuous so the theory of unconditional stability for 
% CN scheme does not apply. You can solve that by reducing theta. The 
% algorithm becomes more stable as theta decreases to 0.