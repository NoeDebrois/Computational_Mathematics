run('A_load_data.m');

%% Robust Frontier 1: 
N = 100;
P1 = Portfolio('AssetList', assetNames);
P1 = setDefaultConstraints(P1); % standard constraints: sum of weights = 1, no short-selling.

P1_RiskSim = zeros(100, N);
P1_RetnSim = zeros(100, N);
P1_Weights = zeros(N, nAssets, 100); %n. sim, n.assets, n portf in a single frontier

for n = 1:N
    R = mvnrnd(ExpRet, CovMatrix);
    NewExpRet = R;
    NewCov = iwishrnd(CovMatrix, nAssets);
    Psim = setAssetMoments(P1, NewExpRet, NewCov);
    w_sim = estimateFrontier(Psim, 100);
    [pf_RiskSim, pf_RetnSim] = estimatePortMoments(Psim, w_sim);
    P1_RiskSim(:, n) = pf_RiskSim;
    P1_RetnSim(:, n) = pf_RetnSim;
    P1_Weights(n, :, :) = w_sim;
end

P1_Weights_mean = squeeze(mean(P1_Weights, 1));

frontier1.id = 3;
frontier1.risk = mean(P1_RiskSim, 2);
frontier1.retn = mean(P1_RetnSim, 2);

% Minimum Variance Portfolio (PortfolioE)
[~, idx_P1_mvp] = min(frontier1.risk);
portfolioE_weights = P1_Weights_mean(:, idx_P1_mvp);
portfolioE = SavedPortfolio("PortfolioE", assetNames, portfolioE_weights, frontier1);
portfolioE.evaluate(Ret);

% Maximum Sharpe Ratio Portfolio (PortfolioG)
P1_sharpe_ratios = frontier1.retn ./ sqrt(frontier1.risk);
[~, idx_P1_msrp] = max(P1_sharpe_ratios);
portfolioG_weights = P1_Weights_mean(:, idx_P1_msrp);
portfolioG = SavedPortfolio("PortfolioG", assetNames, portfolioG_weights, frontier1);
portfolioG.evaluate(Ret);

%% Robust Frontier 2: 
N = 100;
P2 = Portfolio('AssetList', assetNames);
P2 = setDefaultConstraints(P2); % standard constraints: sum of weights = 1, no short-selling.

% Constraint1: Total exposure to sensible sectors has to be less than 50%
P2 = setInequality(P2, sensible_indices, 0.5);

% Constraint2: Total exposure to defensive sectors has to be greater than 30%
P2 = addInequality(P2, -defensive_indices, -0.3);

% Constraint3: Exposure to Communication and Energy sectors has to be between 0.05 and 0.1
lb = zeros(nAssets, 1); % Lower bounds (no short-selling)
ub = ones(nAssets, 1); % Upper bounds (no more than 100% weight)
lb(strcmp(assetNames, 'CommunicationServices')) = 0.05;
ub(strcmp(assetNames, 'CommunicationServices')) = 0.1;
lb(strcmp(assetNames, 'Energy')) = 0.05;
ub(strcmp(assetNames, 'Energy')) = 0.1;
P2 = setBounds(P2, lb, ub);

% Constraint4: Total exposure of cyclical sectors must equal defensive sectors
P2 = setEquality(P2, cyclical_indices-defensive_indices, 0);

P2_RiskSim = zeros(100, N);
P2_RetnSim = zeros(100, N);
P2_Weights = zeros(N, nAssets, 100); %n. sim, n.assets, n portf in a single frontier

for n = 1:N
    R = mvnrnd(ExpRet, CovMatrix);
    NewExpRet = R;
    NewCov = iwishrnd(CovMatrix, nAssets);
    Psim = setAssetMoments(P2, NewExpRet, NewCov);
    w_sim = estimateFrontier(Psim, 100);
    [pf_RiskSim, pf_RetnSim] = estimatePortMoments(Psim, w_sim);
    P2_RiskSim(:, n) = pf_RiskSim;
    P2_RetnSim(:, n) = pf_RetnSim;
    P2_Weights(n, :, :) = w_sim;
end

P2_Weights_mean = squeeze(mean(P2_Weights, 1));

frontier2.id = 4;
frontier2.risk = mean(P2_RiskSim, 2);
frontier2.retn = mean(P2_RetnSim, 2);

% Minimum Variance Portfolio (PortfolioF)
[~, idx_P2_mvp] = min(frontier2.risk);
portfolioF_weights = P2_Weights_mean(:, idx_P2_mvp);
portfolioF = SavedPortfolio("PortfolioF", assetNames, portfolioF_weights, frontier2);
portfolioF.evaluate(Ret);

% Maximum Sharpe Ratio Portfolio (PortfolioH)
P2_sharpe_ratios = frontier2.retn ./ sqrt(frontier2.risk);
[~, idx_P2_msrp] = max(P2_sharpe_ratios);
portfolioH_weights = P2_Weights_mean(:, idx_P2_msrp);
portfolioH = SavedPortfolio("PortfolioH", assetNames, portfolioH_weights, frontier2);
portfolioH.evaluate(Ret);

%% Save
if isfile('portfolios.mat')
    save('portfolios.mat', 'portfolioE', 'portfolioF', 'portfolioG', 'portfolioH', '-append');
else
    save('portfolios.mat', 'portfolioE', 'portfolioF', 'portfolioG', 'portfolioH');
end
disp("A_Q3.m > Portfolios E, F, G, and H saved.");
disp(" ");

%% Plot allocation :
% List of portfolios to visualize
portfolios = {portfolioE, portfolioF, portfolioG, portfolioH};

% Plot asset weights in each portfolio
PortfolioViz.plotWeights(portfolios);

%% Check of the sum of weights :
disp('Portfolio Weight Sums:');
for i = 1:length(portfolios)
    portfolio = portfolios{i};
    total_weight = sum(portfolio.weights);
    fprintf('%s: Total Weight = %.6f\n', portfolio.name, total_weight);
end