function [S,SAV,ST,STAV] = KOU_ASSET_AV(S0, r, T, sigma, p, lambdap, lambdam, lambdaK, M, Nsim)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   cf Simulate_Jump_Diffusion.pdf.
    %   INPUT: 
    %   - S0 = underlying  spot price
    %   - r = risk free rate
    %   - T = time to maturity
    %   - sigma = BM vol 
    %   - p = positive jump prob
    %   - lambdap = positive jump intensity
    %   - lambdam = negative jump intensity
    %   - lambdaK = jump time intensity
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
    Z = randn(Nsim, M+1);
    
    %% Compute drift in Q-dynamics
    % Characteristic exponent for VG (cf VG_NIG.pdf):
    char_exp = @(u) - sigma^2/2 * u.^2 + 1i * u .* lambdaK .* (p ./ (lambdap - 1i * u) - (1 - p) ./ (lambdam + 1i * u)); % Without drift
    % Drift under Q:
    drift = r - char_exp(-1i);

    %% Simulation
    NT = poissrnd(lambdaK * T, Nsim)'; % Column vector of NTs for each sim
    figure;
    for i=1:Nsim
        JumpTimes=sort(rand(NT(i),1) * T);
        for j=1:M
            % 1) Simulate the continuous part of the path:
            X(i, j + 1) = X(i, j) + drift * dt + sigma * sqrt(dt) * Z(i, j);
            XAV(i, j + 1) = XAV(i, j) + drift * dt - sigma * sqrt(dt) * Z(i, j);

            % 2) Add jumps : we use CONDITIONAL SIMULATION of jump times
            % - Loop on all the jump times ;
            % - Check if there are jumps in ((i-1) * dt, i * dt] ;
            % - If yes, add the simulated jump. -> Kou simulated jump.
            for jj=1:NT(i)
                if (JumpTimes(jj)>(j-1)*dt).*(JumpTimes(jj)<=j*dt)
                    u = rand; % RV ~ U((0,1)) to know which type of jump (+ or -).
                    if u < p % Positive jump
                        jumpSize = exprnd(1/lambdap); % jumpsize ~ Exp(lambda_+)
                    else % Negative jump
                        jumpSize = -exprnd(1/lambdam); % jumpsize ~ Exp(lambda_-)
                    end
                    X(i, j + 1) = X(i, j + 1) + jumpSize;
                    XAV(i, j + 1) = XAV(i, j + 1) + jumpSize;
                end
            end
        end
        plot(linspace(0, T, M+1), X(i,:));
        hold on;
        plot(linspace(0, T, M+1), XAV(i,:));
        hold on;
    end
    
    %% Get the paths and the final values
    S = S0 * exp(X);
    SAV = S0 * exp(XAV);
    ST = S(:, end);
    STAV = SAV(:, end);

    %% Check of risk-neutrality
    disp('Risk-Neutrality check. It should be zero:');
    [check, ~, CI_check] = normfit(ST - S0 * exp(r * T));
end