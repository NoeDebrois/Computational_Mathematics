classdef PortfolioViz
    methods (Static)
        function plotWeights(portfolios)
            % Create a new figure
            figure;
        
            % Retrieve asset names from the first portfolio for the legend
            asset_names = portfolios{1}.assetNames;
        
            num_assets = length(asset_names);
        
            % Generate a custom colormap with unique colors for each asset
            cmap = turbo(num_assets);
        
            % Initialize portfolio_names for consistency
            portfolio_names = strings(1, length(portfolios));
            
            if length(portfolios) < 2
                % For a single portfolio, just use its weights
                all_weights = portfolios{1}.weights;
                b = bar(portfolios{1}.name, all_weights, 'stacked');
                portfolio_names(1) = portfolios{1}.name;  % Assign name for single portfolio
            else
                % For multiple portfolios, gather the weights into all_weights
                all_weights = zeros(num_assets, length(portfolios));
                
                % Loop through portfolios to concatenate weights and collect names
                for i = 1:length(portfolios)
                    all_weights(:, i) = portfolios{i}.weights;
                    portfolio_names(i) = portfolios{i}.name;
                end
        
                % Plot the weights as a stacked bar chart
                b = bar(all_weights', 'stacked');
        
                % Set x-ticks and x-tick labels as portfolio names
                xticks(1:length(portfolio_names));
                xticklabels(portfolio_names);
            end
        
            ylabel('Weight');
            title('Asset Weights in Each Portfolio', 'FontSize', 16);
            grid on;
        
            % Set each asset (bar segment) to a different color
            for k = 1:num_assets
                b(k).FaceColor = cmap(k, :);
            end
        
            % Add legend for asset names
            legend(asset_names, 'Location', 'bestoutside', 'FontSize', 16);
        
            % Display the percentages on top of each bar segment
            hold on;
            for i = 1:length(portfolio_names)
                cumulative_weight = 0;
                for j = 1:num_assets
                    % Get the percentage value of the segment
                    percentage = all_weights(j, i) * 100;
                    % Display the percentage on the bar
                    text(i, cumulative_weight + all_weights(j, i) / 2, sprintf('%.1f%%', percentage), ...
                        'Color', 'white', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14, 'FontWeight', 'bold');
                    cumulative_weight = cumulative_weight + all_weights(j, i);
                end
            end
            hold off;
        end

        function plotEfficientFrontier_(portfolios)
            hold on;

            % Initialize a map to store unique frontiers by their id
            uniqueFrontiers = [];
        
            % Loop through each portfolio to plot its efficient frontier
            for i = 1:length(portfolios)
                if ~ isempty(portfolios{i}.frontier)
                    % Extract the current portfolio's frontier
                    frontier = portfolios{i}.frontier;
                    
                    % Check if this frontier has already been plotted by its id
                    frontierID = num2str(frontier.id);
                    if ~ ismember(frontierID, uniqueFrontiers)
                        uniqueFrontiers = [uniqueFrontiers, frontierID];
                        plot(frontier.risk, frontier.retn, 'LineWidth', 1.5, 'DisplayName', strcat("Frontier", frontierID));
                    end
                end
            end
        
            % Set labels and title
            xlabel('Risk (Standard Deviation)');
            ylabel('Expected Return');
            title('Efficient Frontier of Portfolios');
            grid on;

            % Display the legend
            legend('show', 'Location', 'best');

            hold off;
        end

        function plotReturnRisk(portfolios, plot_frontier)
            if nargin < 2
                plot_frontier = false;
            end

            % Create a new figure
            figure;

            % Initialize arrays for risk, return, and portfolio names
            risks = [];
            returns = [];
            portfolio_names = {};
        
            % Extract data from each portfolio
            for i = 1:length(portfolios)
                risks(i) = portfolios{i}.risk;              % Extract risk (x-axis)
                returns(i) = portfolios{i}.retn;  % Extract expected return (y-axis)
                portfolio_names{i} = portfolios{i}.name;    % Extract portfolio name
            end
        
            % Plot return vs risk
            scatter(risks, returns, 100, 'filled', 'DisplayName', 'Portfolios');  % Scatter plot with size 100 for points
        
            % Add labels to each point
            for i = 1:length(portfolios)
                text(risks(i), returns(i), portfolio_names{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
            end

            if plot_frontier
                PortfolioViz.plotEfficientFrontier_(portfolios)
            end

            xlabel('Risk (Standard Deviation)');
            ylabel('Expected Return');
            title('Risk vs Expected Return of Portfolios');
            grid on;
        end

        function plotEfficientFrontier(portfolios)
            % Create a new figure
            figure;
            PortfolioViz.plotEfficientFrontier_(portfolios)
        end

        function plotEvolution(portfolios, dates, Ret)
            figure();

            portfolio_names = strings(1, length(portfolios));
            for i = 1:length(portfolios)
                equity = cumprod(Ret*portfolios{i}.weights);
                equity = 100.*equity/equity(1);
                
                plot(dates(2:end, 1), equity, 'LineWidth', 2)
                hold on;

                portfolio_names(i) = portfolios{i}.name;
            end

            legend(portfolio_names);
        end
    end
end
