run('A_load_data.m');

%% Portfolio object
P = Portfolio('AssetList', assetNames);
P = setDefaultConstraints(P); % standard constraints: sum of weights = 1, no short-selling.
P = estimateAssetMoments(P, LogRet, 'missingdata', false);

%% Compute efficient frontier
pwgt = estimateFrontier(P, 100);
frontier.id = 1;
[frontier.risk, frontier.retn] = estimatePortMoments(P, pwgt);

%% Minimum Variance Portfolio (PortfolioA)
portfolioA_weights = estimateFrontierByRisk(P, min(frontier.risk));
portfolioA = SavedPortfolio("PortfolioA", assetNames, portfolioA_weights, frontier);
portfolioA.evaluate(Ret);

%% Maximum Sharpe Ratio Portfolio (PortfolioB)
portfolioB_weights = estimateMaxSharpeRatio(P);
portfolioB = SavedPortfolio("PortfolioB", assetNames, portfolioB_weights, frontier);
portfolioB.evaluate(Ret);

%% Save
if isfile('portfolios.mat')
    save('portfolios.mat', 'portfolioA', 'portfolioB', '-append');
else
    save('portfolios.mat', 'portfolioA', 'portfolioB');
end
disp("A_Q1.m > Portfolios A and B saved.");
disp(" ");

%% Plot allocation :
% List of portfolios to visualize
portfolios = {portfolioA, portfolioB};

% Plot asset weights in each portfolio
PortfolioViz.plotWeights(portfolios);