classdef SavedPortfolio < handle
    properties
        name
        assetNames
        weights
        frontier
        retn
        risk
        metrics
    end
    
    methods
        function obj = SavedPortfolio(name, assetNames, weights, frontier)
            obj.name = name;
            obj.assetNames = assetNames;
            obj.weights = weights;

            if nargin < 4 || isempty(frontier)
                obj.frontier = [];
            else
                obj.frontier = frontier;
            end
        end

        function ret = evaluate(obj, Ret, varargin)
            % Calculate return and risk
            LogRet = log(Ret);
            ExpRet = mean(LogRet);
            CovMatrix = cov(LogRet);
        
            retn = sum(obj.weights' .* ExpRet);
            risk = sqrt(obj.weights' * CovMatrix * obj.weights);
        
            % Calculate additional metrics
            metrics = obj.getMetrics(Ret);

            ret = {retn, risk, metrics};

            % Parse optional arguments (saveResults and printTable booleans)
            p = inputParser;
            addParameter(p, 'saveResults', true);
            addParameter(p, 'printTable', true);
            parse(p, varargin{:});
        
            % Save results if asked
            if p.Results.saveResults
                obj.retn = retn;
                obj.risk = risk;
                obj.metrics = metrics;
            end
        
            % Print table if asked
            if p.Results.printTable
                disp(obj.name + " - " + "Metrics:");
                disp("Return: " + retn);
                disp("Vol: " + risk);
                disp(" ");

                keys = ["AnnRet", "AnnVol", "Sharpe", "MaxDD", "Calmar", "Herfindahl", "Shannon", "DiversificationRatio", "EntropyVolatility", "EntropyRiskContribution"];
                values = struct2cell(metrics);
                T = array2table(horzcat(values{:}), 'VariableNames', keys);

                disp(T)
                disp(" ");
            end
        end
    end

    methods (Access = private)
        % Private method to calculate all metrics
        function metrics = getMetrics(obj, Ret)
            LogRet = log(Ret);
            CovMatrix = cov(LogRet);

            CumRet = cumprod(Ret*obj.weights);
            CumRet = 100.*CumRet/CumRet(1);

            metrics.annualRtn = obj.annualReturn(CumRet);
            metrics.annualVol = obj.annualVolatility(Ret*obj.weights);
            metrics.sharpe = obj.sharpeRatio(Ret);
            metrics.maxDD = obj.maxDrawdown(CumRet);
            metrics.calmar = obj.calmarRatio(CumRet);
            metrics.h = obj.herfindahlIndex();
            metrics.shannon = obj.shannonIndex();
            metrics.divRatio = obj.diversificationRatio(CovMatrix);
            metrics.entropyVolatility_ = obj.entropyVolatility(CovMatrix);
            metrics.entropyRiskContribution_ = obj.entropyRiskContribution(CovMatrix);
        end

        % Private method to calculate the Annualized Return (CAGR)
        function annualRtn = annualReturn(~, CumRet)
            t = length(CumRet)/252;
            annualRtn = (CumRet(end)/CumRet(1))^(1/t) - 1;
        end
        
        % Private method to calculate the Annualized Volatility
        function annualVol = annualVolatility(~, Ret)
            annualVol = std(Ret) * sqrt(252);
        end
        
        % Private method to calculate the Sharpe Ratio
        function sharpe = sharpeRatio(obj, Ret)
            CumRet = cumprod(Ret*obj.weights);
            CumRet = 100.*CumRet/CumRet(1);

            annualRtn = obj.annualReturn(CumRet);
            annualVol = obj.annualVolatility(Ret*obj.weights);
            sharpe = annualRtn / annualVol;
        end
        
        % Private method to calculate the Maximum Drawdown
        function maxDD = maxDrawdown(~, CumRet)
            maxDD = inf; % Initialize max drawdown
            for i = 2:length(CumRet)
                maxDD = min(maxDD, CumRet(i)/max(CumRet(1:i-1))-1);
            end
        end
        
        % Private method to calculate the Calmar index
        function calmar = calmarRatio(obj, CumRet)
            drawdown = maxDrawdown(obj, CumRet);
            annualRetn = annualReturn(obj, CumRet);
            calmar = annualRetn / abs(drawdown);
        end

        % Private method to calculate the Herfindahl index
        function h = herfindahlIndex(obj)
            h = sum(obj.weights.^2);
        end

        % Private method to calculate Shannon index
        function sIndex = shannonIndex(obj)
            w = obj.weights;
            w = w(w > 0);
            sIndex = -sum(w.*log(w));
        end

        % Private method to calculate the diversification ratio
        function divRatio = diversificationRatio(obj, CovMatrix)
            risk = sqrt(obj.weights' * CovMatrix * obj.weights);
            weighted_volatilities = obj.weights .* sqrt(diag(CovMatrix));
            divRatio = sum(weighted_volatilities) / risk;
        end

        % Private method to calculate entropy based on volatility
        function entropyVolatility = entropyVolatility(obj, CovMatrix)
            volatilities = sqrt(diag(CovMatrix));
            p_vol = ((obj.weights.^2).*(volatilities.^2))/sum((obj.weights.^2).*(volatilities.^2));
            p_vol = p_vol(p_vol > 0);
            entropyVolatility = -sum(p_vol .* log(p_vol));
        end

        % Private method to calculate entropy of risk contribution
        function entropyRiskContribution = entropyRiskContribution(obj, CovMatrix)
            portfolio_vol = sqrt(obj.weights' * CovMatrix * obj.weights);
            marginal_contributions = CovMatrix*obj.weights/portfolio_vol;
            RC = (obj.weights .* marginal_contributions)/portfolio_vol;
            RC = RC(RC > 0);
            entropyRiskContribution = -sum(RC .* log(RC));
        end
    end
end