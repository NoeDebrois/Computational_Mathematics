%% LECTURE 10 - MC for American Options - Noé Debrois - 28/10/2024
% Implementation of Longstaff & Schwartz algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
% Longstaff & Schwartz algorithm --------
% == L=3 ================================
% Quadratic Polynomial ------------------
%
%% === INPUT ============================
S0 = 1;            % Spot price
T = 3;             % Maturity in year(s) 
r = 0.06;          % Risk-free rate
K = 1.1;           % Option strike
M = 3; dt = T / M; % Yearly monitoring
N = 1000;          % Number of simulations
%
%% === SIMULATION OF THE ASSET PRICE ====
% Either we simulate the asset price ourselves or we use the example from
% the Longstaff & Schwartz article.

% -- 1st option : from L&S article :
% S = [1.09 1.08 1.34;
%      1.16 1.26 1.54;
%      1.22 1.07 1.03;
%      0.93 0.97 0.92;
%      1.11 1.56 1.52;
%      0.76 0.77 0.90;
%      0.92 0.84 1.01;
%      0.88 1.22 1.34];

% -- 2nd option : MC simulation : 
sigma = 0.4;
S = AssetBS(r, sigma, S0, T, M, N); % Our asset price model : it can be
% anything : B&S, Kou, etc.
S = S(:, 2:end); % We do not consider t = 0 for early exercise.
%
%% == INITIALIZE =========================
% 1st step : for each option, we set conventionally to T the exercise time
% if the option will never be exercised in the future :
Exercise_Time = M * ones(N, 1);

% Payoff initialization (for the moment, at maturity) :
Put = max(0, K - S(:, end));
%
%% == BACKWARD-IN-TIME ===================
% WARNING :
% We consider only in the money case (Inmoney), i.e K > S_i_j. 
% Because clearly, if K < S_i_j, the exercise would result to a 0 payoff.
% Hence, one doesn't exercise when K < S_i_j.

% Backward-in-time procedure :
for j=M-1:-1:1
    % AT TIME j :
    Inmoney = find(S(:,j) < K); % Get all the index from In The Money Opt.
    S_I = S(Inmoney, j);        % Get all the paths from In The Money Opt.

    % -- Intrinsic Value computation :
    IV = K - S_I; % Payoff if we exercise now (only for In The Money Opt.)

    % -- Continuation Value computation :
    % - Regression : --------------------
    % Basis functions : [1, S, S^2].
    A = [ones(length(S_I), 1), S_I, S_I.^2]; % Here, degree 2 polynomial.

    % b = (b_i)_i = (exp[- r*dt * (j_ex-j)] *  max(K -  S_i_{j_ex},0))_i :
    b = Put(Inmoney) .* exp(-r * dt * (Exercise_Time(Inmoney) - j));

    % Solve A * alpha = b where : 
    % - alpha contains the weights alpha_k_j, k = 1, ..., L=3, L being the
    % degree of the polynomial ;
    % - b = (b_i)_i (cf above) ;
    % - A is rectangular and contains :
    %    - 1 column of ones (of length "length(S_I)") ;
    %    - 1 column = S_I ;
    %    - 1 column = S_I^2.
    alpha = A \ b; 

    % - Continuation Value : ------------
    % Approximation of CV_i_j "=~" SUM_{k=1}^{L} alpha_k_j * S_j_i^{k-1}
    CV = A * alpha; % A * alpha = SUM_{k=1}^{L} A_i_k * alpha_k_j

    % --------

    % IS j AN EXERCISE INSTANT ? 
    % If yes, we need to update the exercise instant for the 
    % corresponding paths.

    % Paths with early exercise at time step j :
    Index = find(IV > CV); % Paths idx where it's better to exercise @ j.
    Early_Exercise = Inmoney(Index); % Paths where it's better to 
    % exercise at time j.

    % Update the price from those options :
    Put(Early_Exercise) = IV(Index); % Their price become their IV.

    % Update the exercise time from those options :
    Exercise_Time(Early_Exercise) = j; % i.e In The Money Options where 
    % it is better to exercise now : we update their exercise time.
end

% Price = Expected value of the discounted payoff : 
[price, ~, CI] = normfit(Put .* exp(-r * dt * Exercise_Time)) 
% [~, EuPrice] = blsprice(S0, K, r, T, sigma)