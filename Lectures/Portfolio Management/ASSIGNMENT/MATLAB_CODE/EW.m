run('A_load_data.m');

%% Portfolio object
P = Portfolio('AssetList', assetNames);
P = setDefaultConstraints(P);
P = estimateAssetMoments(P, LogRet, 'missingdata', false);

%% Compute efficient frontier
pwgt = estimateFrontier(P, 100);
frontier.id = 1;
[frontier.risk, frontier.retn] = estimatePortMoments(P, pwgt);

%% Equally Weighted Portfolio (PortfolioEW)
numAssets = length(assetNames);
equalWeights = ones(numAssets, 1) / numAssets; % Equal weights for all assets
portfolioEW = SavedPortfolio("PortfolioEW", assetNames, equalWeights, frontier);
portfolioEW.evaluate(Ret);
