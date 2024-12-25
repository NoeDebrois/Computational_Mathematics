function [S,SAV,ST,STAV] = MERTON_ASSET_AV(S0, r, T, sigma, lambda, muJ, deltaJ, M, Nsim)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   cf Algorithm 6.2 @ Simulate_Jump_Diffusion.pdf.
    %   INPUT: 
    %   - S0 = underlying  spot price
    %   - r = risk free rate
    %   - T = time to maturity
    %   - sigma = conditional vol 
    %   - lambda = Poisson intensity/rate (jumptimes)
    %   - muJ = Jumpsize ~ Normal(muJ, deltaJ^2)
    %   - deltaJ = Jumpsize ~ Normal(muJ, deltaJ^2)
    %   - M = number of time steps in [0, T]
    %   - Nsim = number of simulations
    %   OUTPUT:
    %   - S = matrix of Nsim underlying paths (one for each row)
    %   - ST = vector of Nsim underlying sim at time T 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % HERE WE CHOOSE TO ONLY NEGATIVELY CORRELATE THE BM (C0 PART) AND USE
    % THE SAME JUMPS. cf page 83 of my written notes for the solution with
    % also the negative correlation of the jumpsizes.

    dt = T / M;
    X = zeros(Nsim, M+1);
    XAV = zeros(Nsim, M+1);
    Z = randn(Nsim, M);
    
    %% Compute drift in Q-dynamics
    % Characteristic exponent for Merton (cf Kou_Merton.pdf):
    char_exp = @(u) - sigma^2/2 * u.^2 + lambda .* (exp(-deltaJ^2/2 * u.^2 + muJ * 1i * u) - 1);
    % Drift under Q:
    drift = r - char_exp(-1i);
    
    %% Simulation
    % - Generates random nb from Poisson distribution of parameter lambda*t:
    NT = poissrnd(lambda * T, Nsim)'; % Column Vector (thx to ')

    for i=1:Nsim
        % - Generates and orders jump times chronologically :
        jumpT = sort(rand(NT(i), 1) * T); % Rescale to [O, T]
        % - Simulate jump sizes, jumpSize ~ Normal(muJ, deltaJ^2) :
        jumpSize = muJ + deltaJ * randn(NT(i), 1);
        for j=1:M
            % 1) Simulate the continuous part of the path :
            X(i, j+1) = X(i, j) + drift * dt + sigma * sqrt(dt) * Z(i, j);
            XAV(i, j+1) = XAV(i, j) + drift * dt - sigma * sqrt(dt) * Z(i, j);

            % 2) Add jumps : we use CONDITIONAL SIMULATION of jump times (cf above)
            % - Loop on all the jump times ;
            % - And check if there are jumps in ((i-1) * dt, i * dt] ;
            % - If yes, add the simulated jump. -> Merton simulated jump.
            for jj=1:NT(i) % Loop on all the jump times
                % Does the jump occur between the 2 timesteps we consider ?
                if jumpT(jj) > (j-1) * dt && jumpT(jj) <= j * dt
                    % If yes, we add the simulated jump size :
                    X(i, j + 1) = X(i, j + 1) + jumpSize(jj);
                    XAV(i, j + 1) = XAV(i, j + 1) + jumpSize(jj);
                end
            end
        end
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