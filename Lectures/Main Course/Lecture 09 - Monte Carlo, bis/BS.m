%% LECTURE 9 - MC Simulation, bis - Noé Debrois - 27/10/2024
% See the explanations at the end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear; close all;
%
%% Geometric Brownian Motion :
% Parameters for geometric Brownian motion
mu = 0.03;          % Drift term
sigma = 0.5;        % Volatility term
S0 = 1;             % Initial asset price
T = 10;             % Time horizon
disp('Computes the exact expected value of S(T) using the formula: E[S(T)] = S_0 * exp(mu*T) :');
% Calculate and display the exact expected value of S(T)
Exact = S0 * exp(mu * T)
disp('****----****----****----****')
%
%% Monte Carlo simulation :
disp('Exact Solution MC Simulation :')
Nsim = 1e6; % Number of MC simulations
% Simulate S(T) using exact formula and calculate mean and 95% CI
[Value, ~, CI] = normfit( ...
    S0 * exp((mu - sigma^2 / 2) * T + sigma * sqrt(T) * randn(Nsim, 1)))
disp('****----****----****----****')
%
%% Euler scheme simulations :
% Euler scheme with 1 time step (M=1)
M = 1;              % Number of time steps
dt = T / M;         % Time increment
S = S0;             % Initialize asset price at time 0
disp('Euler scheme with 1 time step :')
% Simulate asset path using Euler scheme with M=1
for i = 1:M
    S = S .* (1 + mu * dt + sigma * sqrt(dt) * randn(Nsim, 1));
end
% Calculate mean and 95% CI of simulated S(T) using normfit
[Value, ~, CI] = normfit(S)
disp('****----****----****----****')

% Euler scheme with 10 time steps (M=10)
M = 10;             % Number of time steps
dt = T / M;         % Time increment
S = S0;             % Reset asset price for new simulation
disp('Euler scheme with 10 time steps :')
% Simulate asset path using Euler scheme with M=10
for i = 1:M
    S = S .* (1 + mu * dt + sigma * sqrt(dt) * randn(Nsim, 1));
end
% Calculate mean and 95% CI of simulated S(T) using normfit
[Value, ~, CI] = normfit(S)
disp('****----****----****----****')

% Euler scheme with 100 time steps (M=100)
M = 100;            % Number of time steps
dt = T / M;         % Time increment
S = S0;             % Reset asset price for new simulation
disp('Euler scheme with 100 time steps :')
% Simulate asset path using Euler scheme with M=100
for i = 1:M
    S = S .* (1 + mu * dt + sigma * sqrt(dt) * randn(Nsim, 1));
end
% Calculate mean and 95% CI of simulated S(T) using normfit
[Value, ~, CI] = normfit(S)
%
%% Explanations :
% This MATLAB code simulates the expected value of an asset price S(T)
% at a future time T using a geometric Brownian motion model, with a 
% given drift mu and volatility sigma. The code implements both the exact 
% solution for E[S(T)] and simulations using the Euler scheme with varying
% time steps to approximate the solution.
% 
% 1. Initialization and Exact Solution Calculation:
%    - Sets parameters for drift (mu=0.03), volatility (sigma=0.5), 
%      initial price (S0=1), and time horizon (T=10).
%    - Computes the exact expected value of S(T) using the formula:
%      E[S(T)] = S_0 * exp(mu*T)
% 
% 2. Exact Solution Simulation Using Monte Carlo:
%    - Simulates N_{sim} = 10^6 asset paths using the exact solution:
%      S(T) = S_0 * exp([mu - sigma^2/2] * T + sigma * sqrt{T} * Z}
%      where Z is a standard normal random variable ('randn(Nsim,1)').
%    - Uses 'normfit' to estimate the mean and 95% confidence interval (CI)
%      of the simulated S(T) values.
% 
% 3. Euler Scheme Approximation:
%    - Implements the Euler-Maruyama scheme to approximate S(T) with three 
%      different time steps (M=1, M=10, M=100), where M is the number of 
%      subintervals within T.
%    - For each value of M :
%      - Computes the time increment (dt = T / M).
%      - Initializes S to the starting value (S0) and iteratively updates S
%        according to:
%        S = S * (1 + mu * dt + sigma * sqrt{dt} * Z)
%      - Uses 'normfit' to compute the mean and 95% CI of the simulated 
%        S(T) values for each time step. 
% 
% This code essentially compares the exact expected value of S(T) with 
% simulated values using the Euler scheme, observing the effect of 
% different time steps on the accuracy of the simulation.