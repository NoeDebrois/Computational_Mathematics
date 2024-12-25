run("A_load_data.m");

%% Building the views
v = 2;
tau = 1/nAssets;
P = zeros(v, nAssets);
q = zeros(v, 1);
Omega = zeros(v);

% View 1: Energy annual performance of 3%
P(1, strcmp(assetNames, 'Energy')) = 1;
q(1) = 0.03;

% View 2: Quality outperforms Momentum by 1%
P(2, strcmp(assetNames, 'Quality')) = 1;
P(2, strcmp(assetNames, 'Momentum')) = -1;
q(2) = 0.01;

Omega(1,1) = tau.*P(1,:)*CovMatrix*P(1,:)';
Omega(2,2) = tau.*P(2,:)*CovMatrix*P(2,:)';

% from annual view to daily view
bizyear2bizday = 1/252;
q = q*bizyear2bizday;
Omega = Omega*bizyear2bizday;

%% market implied ret
wCap = caps/sum(caps);
lambda = 1.2;
mu_mkt = lambda .* CovMatrix * wCap';
C = tau .* CovMatrix;

%% Black Litterman
muBL = inv(inv(C)+P'*inv(Omega)*P)*(P'*inv(Omega)*q + inv(C)*mu_mkt); 
covBL = inv(P'*inv(Omega)*P + inv(C));

P = Portfolio('NumAssets', nAssets, 'Name', 'MV with BL');
P = setDefaultConstraints(P);
P = setAssetMoments(P, muBL, CovMatrix+covBL);

%% Compute efficient frontier
pwgt = estimateFrontier(P, 100);
frontier.id = 5;
[frontier.risk, frontier.retn] = estimatePortMoments(P, pwgt);

%% Minimum Variance Portfolio (PortfolioI)
portfolioI_weights = estimateFrontierByRisk(P, min(frontier.risk));
portfolioI = SavedPortfolio("PortfolioI", assetNames, portfolioI_weights, frontier);
portfolioI.evaluate(Ret);

%% Maximum Sharpe Ratio Portfolio (PortfolioL)
portfolioL_weights = estimateMaxSharpeRatio(P);
portfolioL = SavedPortfolio("PortfolioL", assetNames, portfolioL_weights, frontier);
portfolioL.evaluate(Ret);

%% Save
if isfile('portfolios.mat')
    save('portfolios.mat', 'portfolioI', 'portfolioL', '-append');
else
    save('portfolios.mat', 'portfolioI', 'portfolioL');
end
disp("A_Q4.m > Portfolios I and L saved.");
disp(" ");

%% Plot allocation :
% List of portfolios to visualize
portfolios = {portfolioI, portfolioL};

% Plot asset weights in each portfolio
PortfolioViz.plotWeights(portfolios);