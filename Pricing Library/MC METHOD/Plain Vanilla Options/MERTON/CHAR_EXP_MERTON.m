function V = CHAR_EXP_MERTON(v, parameters)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic exponent of logprices in Merton model
    % INPUT: - parameters = parameters structure{
    %                     - parameters.sigma = BM vol 
    %                     - parameters.mu = jump drift
    %                     - parameters.delta = jump vol
    %                     - parameters.lambdaK = jump time intensity
    %                     }
    % OUTPUT: characteristic exponent for the MERTON model under risk neutral measure 
    % cf KOU_MERTON.pdf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % De-struct parameters 
    sigma = parameters.sigma;
    mu = parameters.mu;
    delta = parameters.delta;
    lambdaK  = parameters.lambdaK;
    
    % Evaluate char exp under risk-neutral measure 
    V = @(v) - sigma^2/2 * v.^2 + lambdaK .* (exp(-delta^2/2 * v.^2 + mu * 1i * v) - 1); % Without drift
    drift_rn = -V(-1i); % Drift Risk_neutral
    V = drift_rn * 1i * v + V(v);
end