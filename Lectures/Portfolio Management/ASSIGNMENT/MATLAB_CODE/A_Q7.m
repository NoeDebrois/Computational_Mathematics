run('A_load_data.m');

% Function to minimize
function SR_ES = sharpeRatioES(x, ExpRet, LogRet, alpha)
    % Portfolio Expected Return
    portReturn = ExpRet*x;

    % Portfolio Losses
    portLosses = -LogRet*x;

    % Expected Shortfall (ES)
    sortedLosses = sort(portLosses, 'ascend');
    varIndex = floor(alpha * length(sortedLosses)); % Index for VaR
    ES = mean(sortedLosses(varIndex:length(sortedLosses))); % Average of worst losses

    % Sharpe Ratio
    SR_ES = portReturn/ES;

    % Negative for minimization
    SR_ES = -SR_ES;
end

% Standard constraints
x0 = rand(nAssets, 1);
x0 = x0./sum(x0);
lb = zeros(1, nAssets);
ub = ones(1, nAssets);
Aeq = ones(1, nAssets);
beq = 1;

[portfolioQ_weights, ~] = fmincon(@(x) sharpeRatioES(x, ExpRet, LogRet, 0.95), x0, [], [], Aeq , beq, lb, ub);
portfolioQ = SavedPortfolio("PortfolioQ", assetNames, portfolioQ_weights);
portfolioQ.evaluate(Ret);

%% Save
if isfile('portfolios.mat')
    save('portfolios.mat', 'portfolioQ', '-append');
else
    save('portfolios.mat', 'portfolioQ');
end
disp("A_Q7.m > Portfolio Q saved.");
disp(" ");

%% Plot allocation :
% List of portfolios to visualize
portfolios = {portfolioQ};

% Plot asset weights in each portfolio
PortfolioViz.plotWeights(portfolios);