run('A_load_data.m');

%% Portfolio object
P = Portfolio('AssetList', assetNames);
P = setDefaultConstraints(P); % standard constraints: sum of weights = 1, no short-selling.

% Constraint1: Total exposure to sensible sectors has to be less than 50%
P = setInequality(P, sensible_indices, 0.5);

% Constraint2: Total exposure to defensive sectors has to be greater than 30%
P = addInequality(P, -defensive_indices, -0.3);

% Constraint3: Exposure to Communication and Energy sectors has to be between 0.05 and 0.1
lb = zeros(nAssets, 1); % Lower bounds (no short-selling)
ub = ones(nAssets, 1); % Upper bounds (no more than 100% weight)
lb(strcmp(assetNames, 'CommunicationServices')) = 0.05;
ub(strcmp(assetNames, 'CommunicationServices')) = 0.1;
lb(strcmp(assetNames, 'Energy')) = 0.05;
ub(strcmp(assetNames, 'Energy')) = 0.1;
P = setBounds(P, lb, ub);

% Constraint4: Total exposure of cyclical sectors must equal defensive sectors
P = setEquality(P, cyclical_indices-defensive_indices, 0);

P = estimateAssetMoments(P, LogRet,'missingdata',false);

%% Compute efficient frontier
pwgt = estimateFrontier(P, 100);
frontier.id = 2;
[frontier.risk, frontier.retn] = estimatePortMoments(P, pwgt);

%% Minimum Variance Portfolio (PortfolioC)
portfolioC_weights = estimateFrontierByRisk(P, min(frontier.risk));
portfolioC = SavedPortfolio("PortfolioC", assetNames, portfolioC_weights, frontier);
portfolioC.evaluate(Ret);

%% Maximum Sharpe Ratio Portfolio (PortfolioD)
portfolioD_weights = estimateMaxSharpeRatio(P);
portfolioD = SavedPortfolio("PortfolioD", assetNames, portfolioD_weights, frontier);
portfolioD.evaluate(Ret);

%% Save
if isfile('portfolios.mat')
    save('portfolios.mat', 'portfolioC', 'portfolioD', '-append');
else
    save('portfolios.mat', 'portfolioC', 'portfolioD');
end
disp("A_Q2.m > Portfolios C and D saved.");
disp(" ");

%% Plot allocation :
% List of portfolios to visualize
portfolios = {portfolioC, portfolioD};

% Plot asset weights in each portfolio
PortfolioViz.plotWeights(portfolios);