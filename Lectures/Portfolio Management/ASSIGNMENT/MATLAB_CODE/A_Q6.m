run('A_load_data.m');

stdLogRet = std(LogRet);
LogRetStd = (LogRet-ExpRet)./stdLogRet;

%% PCA used to find the minimum number of factors that explains more than the 95% of the cumulative variance
k = nAssets;
[~, ~, ~, ~, explained, ~] = pca(LogRetStd, 'NumComponents', k);

% Calculate cumulative explained variance :
CumExplainedVar = cumsum(explained)/100;

% Minimum number of components explaining more than 95% of the variance :
minFactors = find(CumExplainedVar >= 0.95, 1); % 1st index where the condition is met.

% Print the result
disp(['Number of factors explaining more than 95% of the variance : ', num2str(minFactors)]);
disp(['Cumulative explained variance with ', num2str(minFactors), ' factors : ', num2str(CumExplainedVar(minFactors))]);

% % Plot CumExplainedVar
% f = figure();
% plot(1:nAssets, CumExplainedVar, 'm')
% yline(0.95, '--r', 'LineWidth', 1.2, 'Label', '95% Threshold');
% hold on
% scatter(1:nAssets, CumExplainedVar,'m', 'filled')
% grid on
% title('Total Percentage of Explained Variances for the first n-components')
% xlabel('Total number of Principal Components')
% ylabel('Percentage of Explained Variances')
% legend('Explained Variance', '95% Threshold', 'Location', 'southeast');

%% Real PCA
k = minFactors;
[factorLoading, factorRetn, latent, r, explained, mu] = pca(LogRetStd, 'NumComponents', k);

% Reconstruct asset returns
reconReturnStd = factorRetn*factorLoading' + mu;
reconReturn = stdLogRet.*reconReturnStd + ExpRet; % 'Unstandardize'

% Reconstruct covar matrix
covarFactor = cov(factorRetn);

unexplainedRetn = LogRetStd - reconReturnStd;   % To compute D matrix
unexplainedCovar = diag(cov(unexplainedRetn));
D = diag(unexplainedCovar);

covarAssetStd = factorLoading*covarFactor*factorLoading' + D;
covarAsset = diag(stdLogRet.^2)*covarAssetStd; % 'Unstandardize'

%% Optimization max(ret)
% Function to minimize
func = @(x) -ExpRet*x;

% Standard constraints
x0 = rand(nAssets, 1);
x0 = x0./sum(x0);
lb = zeros(1, nAssets);
ub = ones(1, nAssets);
Aeq = ones(1, nAssets);
beq = 1;

% Stricter options to try to converge :
[portfolioP_weights, ~] = fmincon(func, x0, [], [], Aeq, beq, lb, ub, @(x) volatilityConstraint(x, covarAsset, 0.75));
portfolioP = SavedPortfolio("PortfolioP", assetNames, portfolioP_weights);
portfolioP.evaluate(Ret);

%% Save
if isfile('portfolios.mat')
    save('portfolios.mat', 'portfolioP', '-append');
else
    save('portfolios.mat', 'portfolioP');
end
disp("A_Q6.m > Portfolio P saved.");
disp(" ");

%% Functions
% Volatility constraint
function [c, ceq] = volatilityConstraint(x, covar, volTarget)
    c = sqrt(x'*covar*x) - volTarget;
    ceq = 0;
end