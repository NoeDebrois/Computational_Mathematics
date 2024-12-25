function V = CHAR_EXP_NIG(v, parameters)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic exponent of logprices in VG model
    % INPUT: - parameters = parameters structure{
    %                     - parameters.sigma = volatility of the BM
    %                     - parameters.theta = drift of the BM
    %                     - parameters.k = variance of the subordinator
    %                     }
    % cf VG_NIG.pdf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % De-struct parameters
    sigma = parameters.sigma;
    theta = parameters.theta;
    k     = parameters.kNIG;
    
    % Evaluate char exp under risk-neutral measure
    V = @(v) (1 - sqrt(1 + v.^2 .* sigma^2 * k - 2 * 1i * theta * v * k)) / k; % Without drift
    drift_rn = -V(-1i); % Drift Risk_neutral
    V = drift_rn * 1i * v + V(v);
end