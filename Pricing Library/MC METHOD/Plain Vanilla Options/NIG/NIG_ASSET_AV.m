function [S, SAV, ST, STAV] = NIG_ASSET_AV(S0, r, T, sigma, theta, kNIG, M, Nsim)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   ALGORITHM 6.11 @ Simulate_VG_NIG.pdf.
    %   INPUT: 
    %   - S0 = underlying  spot price
    %   - r = risk free rate
    %   - T = time to maturity
    %   - sigma = conditional vol, volatility of the BM
    %   - theta = conditional part of the drift, drift of the BM
    %   - kNIG  = parameter k, variance of the subordinator
    %   - M = number of time steps in [0, T]
    %   - Nsim = number of simulations
    %   OUTPUT:
    %   - S = matrix of Nsim underlying paths (one for each row)
    %   - ST = vector of Nsim underlying sim at time T 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt = T / M;
    X = zeros(Nsim, M+1);
    XAV = zeros(Nsim, M+1);
    Z = randn(Nsim, M+1);
    
    %% Compute drift in Q-dynamics
    % Characteristic exponent for VG (cf VG_NIG.pdf):
    char_exp = @(u) (1 - sqrt(1 + u.^2 .* sigma^2 * kNIG - 2 * 1i * theta * u * kNIG)) / kNIG; % Without drift
    % Drift under Q:
    drift = r - char_exp(-1i);
    
    %% Simulation
    for i=1:M
        % STEP 1:
        dS = icdf('InverseGaussian', rand(Nsim, 1), dt, dt^2 / kNIG); % Subordinator
        % STEP 2: (already done above) simulate N(0, 1) RV: Z(:, i).
        X(:, i+1) = X(:, i) + drift * dt + sigma * sqrt(dS) .* Z(:, i) +...
            theta * dS;
        XAV(:, i+1) = XAV(:, i) + drift * dt - sigma * sqrt(dS) .* Z(:, i) +...
            theta * dS;
    end
    
    %% Get the paths and the final values
    S = S0 * exp(X);
    SAV = S0 * exp(XAV);
    ST = S(:, end);
    STAV = SAV(:, end);

    %% Check of risk-neutrality
    disp('Risk-Neutrality check. It should be zero:')
    [check, ~, CI_check] = normfit(ST - S0 * exp(r * T))
end