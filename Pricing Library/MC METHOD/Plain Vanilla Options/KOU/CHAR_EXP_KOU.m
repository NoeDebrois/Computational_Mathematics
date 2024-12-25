function V = CHAR_EXP_KOU(v, parameters)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic exponent of logprices in BS model
    % INPUT: - parameters = parameters structure{
    %                     - parameters.sigma = BM vol 
    %                     - parameters.p = prob of positive jump
    %                     - parameters.lambdap = pos jump intensity
    %                     - parameters.lambdam = neg jump intensity
    %                     - parameters.lambdaK = jump time intensity
    %                     }
    % OUTPUT: characteristic exponent for the KOU model under risk neutral measure 
    % cf KOU_MERTON.pdf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % De-struct parameters 
    sigma = parameters.sigma;
    p = parameters.p;
    lambdap = parameters.lambdap;
    lambdam = parameters.lambdam;
    lambdaK  = parameters.lambdaK;
    
    % Evaluate char exp under risk-neutral measure 
    V = @(v) - sigma^2/2 * v.^2 + 1i * v .* lambdaK .* (p ./ (lambdap - 1i * v) - (1 - p) ./ (lambdam + 1i * v)); % Without drift
    drift_rn = -V(-1i); % Drift Risk_neutral
    V = drift_rn * 1i * v + V(v);
end