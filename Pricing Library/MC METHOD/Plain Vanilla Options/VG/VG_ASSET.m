function [S,ST] = VG_ASSET(S0, r, T, sigma, theta, kVG, M, Nsim)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   ALGORITHM 6.11 @ Simulate_VG_NIG.pdf.
    %   INPUT: 
    %   - S0 = underlying  spot price
    %   - r = risk free rate
    %   - T = time to maturity
    %   - sigma = conditional vol 
    %   - theta = conditional part of the drift
    %   - kVG   = parameter k
    %   - M = number of time steps in [0, T]
    %   - Nsim = number of simulations
    %   OUTPUT:
    %   - S = matrix of Nsim underlying paths (one for each row)
    %   - ST = vector of Nsim underlying sim at time T 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt = T / M;
    X = zeros(Nsim, M+1);
    Z = randn(Nsim, M+1);
    
    %% Compute drift in Q-dynamics
    % Characteristic exponent for VG (cf VG_NIG.pdf):
    char_exp = @(u) -log(1 + u.^2 * sigma^2 * kVG / 2 - 1i * theta * kVG * u) / kVG;
    % Drift under Q:
    drift = r - char_exp(-1i);
    
    %% Simulation
    for i=1:M
        % STEP 1:
        dS = kVG * icdf('gamma', rand(Nsim, 1), dt / kVG, 1);
        % STEP 2: (already done above) simulate N(0, 1) RV: Z(:, i).
        X(:, i+1) = X(:, i) + drift * dt + sigma * sqrt(dS) .* Z(:, i) +...
            theta * dS;
    end
    
    %% Get the paths and the final values
    S = S0 * exp(X);
    ST = S(:, end);

    %% Check of risk-neutrality
    disp('Risk-Neutrality check. It should be zero:')
    [check, ~, CI_check] = normfit(ST - S0 * exp(r * T))
end