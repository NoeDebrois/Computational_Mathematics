function [price] = BS(vol, spot, strike, rf, ttm)

%--------------------------------------------------------------------------

% Price EU Call with BS for multiple (K, sigma(K)) and maturity at once

% INPUT
% vol >> volatility
% spot
% strike
% rf >> risk free rate
% maturity

% OUTPUT
% price

%--------------------------------------------------------------------------

price = zeros(size(vol));
for i = 1:length(ttm)
    % Compute N(d1) and N(d2)
    d1  = (log(spot./strike)+(rf(i)+vol(i,:).^2*0.5)*ttm(i))./(vol(i,:)*sqrt(ttm(i)));
    d2  = d1-vol(i,:)*sqrt(ttm(i));
    nd1 = normcdf(d1); 
    nd2 = normcdf(d2);
    
    % Compute price
    price(i,:) = spot*nd1-strike.*exp(-rf(i)*ttm(i)).*nd2;
end

end