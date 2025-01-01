function plot_vol_3d(strikes, ttm, vol_BS, vol_HJM)

%--------------------------------------------------------------------------

% Function to plot 3D surfaces for Black-Scholes and HJM model implied
% volatility

% INPUT
% strikes >> Array of strike prices (1 x m)
% ttm >> Array of time to maturities (1 x n)
% prices_BS >> n x m matrix of Black-Scholes prices
% prices_HJM >> n x m matrix of HJM model prices

%--------------------------------------------------------------------------

% Create a grid for strikes and TTM
[K_grid, TTM_grid] = meshgrid(strikes, ttm);

% Create the figure
figure();

% Plot Black-Scholes prices
surf(K_grid, TTM_grid, vol_BS, 'FaceColor', 'blue', 'FaceAlpha', 0.6, 'EdgeColor', 'black');
hold on;

% Plot HJM model prices
surf(K_grid, TTM_grid, vol_HJM, 'FaceColor', 'red', 'FaceAlpha', 0.6, 'EdgeColor', 'black');

% Add labels, legend, and title
xlabel('Strike');
ylabel('TTM');
zlabel('Implied volatility');
legend('Black-Scholes', 'HJM');
grid on;

% Set a view angle for better visualization
view(45, 30);

hold off;

end
