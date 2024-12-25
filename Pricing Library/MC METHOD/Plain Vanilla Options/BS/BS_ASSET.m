function [S,ST] = BS_ASSET(S0, r, T, sigma, M, Nsim)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   cf Simulate_Jump_Diffusion.pdf.
    %   INPUT: 
    %   - S0 = underlying  spot price
    %   - r = risk free rate
    %   - T = time to maturity
    %   - sigma = BM vol 
    %   - M = number of time steps in [0, T]
    %   - Nsim = number of simulations
    %   OUTPUT:
    %   - S = matrix of Nsim underlying paths (one for each row)
    %   - ST = vector of Nsim underlying sim at time T 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt = T / M;
    
    %% Compute drift in Q-dynamics
    % Drift under Q:
    drift = r - sigma^2 / 2;

    %% Simulation
    X = zeros(Nsim, M + 1); 
    Z = randn(Nsim, M);

    for i=1:Nsim
        for j=1:M
            % 1) Simulate the continuous part of the path:
            X(i, j + 1) = X(i, j) + drift * dt + sigma * sqrt(dt) * Z(i, j);
        end
    end
    
    %% Get the paths and the final values
    S = S0 * exp(X);
    ST = S(:, end);

    figure;
    for i=1:Nsim
        plot(linspace(0, T, M+1), S(i,:));
        hold on;
    end    

    %% Check of risk-neutrality
    disp('Risk-Neutrality check. It should be zero:');
    [check, ~, CI_check] = normfit(ST - S0 * exp(r * T))
end
