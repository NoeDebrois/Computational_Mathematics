%% Simulation of the Merton model :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
% S_T = S0 * exp(X_T) with :
% (X_t)_t ~ Merton model : Jumpsize ~ Normal(muJ,deltaJ^2).
%
% General parameters :
S0 = 1; % Initial value
mu = 0.05; % Drift (as in GBM)
sigma = 0.4; % Volatility (as in GBM)
lambda = 2; % Poisson intensity/rate
% Jumpsize parameters : Jumpsize ~ Normal(muJ, deltaJ^2) :
muJ = 0.01; % mean of the size of the jumps
deltaJ = 0.2; % std deviation of the size of the jumps
% Simulation parameters :
T = 2; % Maturity
M = 100; % Number of steps in time
dt = T/M; % Time step
%% Simulation :
% Initialisations of X (in the exp) & Z (vector of N(0, 1) RV) :
X = zeros(M+1,1);
Z = randn(M,1);
%
% CONDITIONAL SIMULATION of jump times (the easiest way) :
% cf "Simulate_Jump_Diffusion.pdf" : ALGORITHM 6.2
% - Generates random nb from Poisson distribution of parameter lambda*t:
NT = poissrnd(lambda*T);
% - Generates and orders jump times chronologically :
jumpT = sort(rand(1,NT)*T); % sort(a vector of U() RV of length NT, rescaled between 0 and T)
% - Simulate jump sizes, jumpSize ~ Normal(muJ, deltaJ^2) :
jumpSize = muJ + deltaJ * randn(NT,1);
% 
for i=1:M
    % 1) Simulate the continuous part of the path :
    %
    X(i+1) = X(i) + mu * dt + sigma * sqrt(dt) * Z(i);
    %
    % 2) Add jumps : we use CONDITIONAL SIMULATION of jump times (cf above)
    % - Loop on all the jump times ;
    % - And check if there are jumps in ((i-1) * dt, i * dt] ;
    % - If yes, add the simulated jump. -> Merton simulated jump.
    for j=1:NT % Loop on all the jump times
        % Does the jump occur between the 2 timesteps we consider ?
        if jumpT(j) > (i-1) * dt && jumpT(j) <= i * dt
            % If yes, we add the simulated jump size :
            X(i+1) = X(i+1) + jumpSize(j);
        end
    end
end
%% Plot :
figure;
plot(S0 * exp(X));
title("One simulated path via Merton model :")