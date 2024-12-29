% Longstaff & Schwartz algorithm, with L = 3 (quadratic polynomials):
% For Bermudan options (i.e American but with discrete monitoring).

%% Parameters:
S0 = 1;            % Spot price
T = 3;             % Maturity in year(s) 
r = 0.06;          % Risk-free rate
sigma = 0.4;       % Volatility
K = 1.1;           % Option strike
M = 3; dt = T / M; % Yearly monitoring
N = 10000;          % Number of simulations

%% Asset price simulation (we can use any Lévy):
[S, ~] = BS_ASSET(S0, r, T, sigma, M, N);
S = S(:, 2:end); % We do not consider t = 0 for early exercise.

%% Initializations:
% Payoff at maturity initialization:
Put = max(K - S(:, end), 0);
% Exercise times initialized at maturity:
Exercise_Time = M * ones(N, 1);

%% Backward-in-time procedure:
for j=M-1:-1:1
    % At time j:
    Inmoney = find(S(:,j) < K); % Get all the index from ITM Opt.
    S_I = S(Inmoney, j);        % Get asset values at j from ITM Opt.
    
    % Intrinsic value for itm options:
    IV = K - S_I; % Payoff if exercise NOW

    % Continuation value for itm options: alpha = A^{-1}*b
    A = [ones(length(S_I), 1), S_I, S_I.^2]; % rows=nb of itm options,
                                             % cols=S_j^k.
    b = Put(Inmoney) .* exp(-r * dt * (Exercise_Time(Inmoney) - j)); % Discounted payoff at exercise time.
    
    % Find the optimized weights at time j:
    alpha = A \ b; 

    % Compute CV_j:
    CV = A * alpha;

    % Is j an exercise instant ? 
    Index = find(IV > CV); % Paths idx where it's better to exercise @ j.
    Early_Exercise = Inmoney(Index); % Index of ITM options w/ early exercise.
    
    % Update exercise time for early exercise ITM options:
    Exercise_Time(Early_Exercise) = j;

    % Update the price of the early exercised options:
    Put(Early_Exercise) = IV(Index);
end

%% Get the final price:
% Price = Expected value of the discounted payoff : 
[price, ~, CI] = normfit(Put .* exp(-r * dt * Exercise_Time))