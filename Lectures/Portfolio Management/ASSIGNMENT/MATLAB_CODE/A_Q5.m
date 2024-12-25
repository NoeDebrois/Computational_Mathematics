run('A_load_data.m');

%% Constraints definition
% Sum of weights = 1
Aeq = ones(1, nAssets);
beq = 1;

% Weights between 0 and 1
lb = zeros(1, nAssets);
ub = ones(1, nAssets);

% Defensive sector counts for less than 20%
A = defensive_indices;
b = 0.2;

%% Capitalization weighted portfolio (PortfolioCap)
portfolioCap_weights = caps'/sum(caps);
portfolioCap = SavedPortfolio("PortfolioCap", assetNames, portfolioCap_weights);
portfolioCap.evaluate(Ret);

%% Most Diversified Portfolio (PortfolioM)
x0 = 1/nAssets*ones(nAssets, 1); % equally weighted portfolio
[portfolioM_weights, ~] = fmincon(@(x) -getDiversificationRatio(x, LogRet), x0, A, b, Aeq, beq, lb, ub, @(x) sumAbsDiffConstraint(x, portfolioCap_weights, 0.3));
portfolioM = SavedPortfolio("PortfolioM", assetNames, portfolioM_weights);
portfolioM.evaluate(Ret);

%% Max Entropy in Risk Contributions (PortfolioN)
x0 = zeros(nAssets ,1);
x0(1, 1) = 1;
[portfolioN_weights, ~] = fmincon(@(x) -getEntropy(getRiskContributions(x, LogRet)), x0, A, b, Aeq, beq, lb, ub, @(x) sumAbsDiffConstraint(x, portfolioCap_weights, 0.3));
portfolioN = SavedPortfolio("PortfolioN", assetNames, portfolioN_weights);
portfolioN.evaluate(Ret);

%% Save
if isfile('portfolios.mat')
    save('portfolios.mat', 'portfolioCap', 'portfolioM', 'portfolioN', '-append');
else
    save('portfolios.mat', 'portfolioCap', 'portfolioM', 'portfolioN');
end
disp("A_Q5.m > Portfolios Cap, M, and N saved.");
disp(" ");

%% Plot allocation :
% List of portfolios to visualize
portfolios = {portfolioCap, portfolioM, portfolioN};

% Plot asset weights in each portfolio
PortfolioViz.plotWeights(portfolios);

%% Functions definition
function [c, ceq] = sumAbsDiffConstraint(x, wCap, diffTolerance)
    c = [];
    ceq = sum(abs(x - wCap)) - diffTolerance;
end

function DR = getDiversificationRatio(x, LogRet)
    % Calculate individual asset volatilities
    sigma_i = std(LogRet);  % 1 x N row vector of asset volatilities

    % Portfolio volatility
    Sigma = cov(LogRet);     % N x N covariance matrix of returns
    portfolio_vol = sqrt(x' * Sigma * x); % Portfolio standard deviation

    % Diversification ratio
    DR = sum(x' .* sigma_i)/portfolio_vol;
end

function H = getEntropy(x)
    % Remove any zero weights to avoid issues with log(0)
    x_nonzero = x(x > 0);
    
    % Calculate entropy
    H = -sum(x_nonzero .* log(x_nonzero));
end

function RC = getRiskContributions(x, LogRet)
    % Covariance matrix of asset returns
    Sigma = cov(LogRet);  % N x N covariance matrix
    
    % Portfolio volatility
    portfolio_vol = sqrt(x' * Sigma * x); % Portfolio standard deviation
    
    % Calculate the marginal risk contributions
    marginal_contributions = Sigma*x/portfolio_vol; % N x 1 vector of marginal contributions
    
    % Calculate risk contributions
    RC = (x .* marginal_contributions)/portfolio_vol; % N x 1 vector of risk contributions
end