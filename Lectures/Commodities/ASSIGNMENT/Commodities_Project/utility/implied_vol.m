function implied_vol = implied_vol(priceMatrix, S, K, r, T)

%--------------------------------------------------------------------------

% Function to compute the implied volatility matrix from a given price matrix
% using the Black-Scholes model.

% Inputs:
% priceMatrix >> n x m matrix of option prices
% S >> Spot price (scalar)
% K >> Strike price (scalar or n x m matrix)
% r >> Risk-free rate (scalar)
% T >> Time to maturity (scalar or n x m matrix)
% optionType >> 'call' or 'put' (string)
%
% Output:
% implied_vol >> n x m matrix of implied volatilities

%--------------------------------------------------------------------------

% Dimensions of the price matrix
[n, m] = size(priceMatrix);

% Initialize the implied volatility matrix
implied_vol = zeros(n, m);

% Nested function to calculate the Black-Scholes price
function price = black_scholes(S, K, r, T, sigma)
    d1 = (log(S ./ K) + (r + 0.5 * sigma.^2) .* T) ./ (sigma .* sqrt(T));
    d2 = d1 - sigma .* sqrt(T);
    price = S .* normcdf(d1) - K .* exp(-r .* T) .* normcdf(d2);
end

% Loop over each element in the price matrix
for i = 1:n
    for j = 1:m
        objFun = @(sigma) black_scholes(S, K(j), r(i), T(i), sigma) - priceMatrix(i,j);
        implied_vol(i,j) = fzero(objFun, 0.5); 
    end
end
end
