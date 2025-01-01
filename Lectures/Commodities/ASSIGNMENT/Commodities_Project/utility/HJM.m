function [price] = HJM(A, n, S0, r, T, strike, params, model, flag)

%--------------------------------------------------------------------------

% Price call options for different strikes and maturities using CM formula
% with FFT and different Levy models

% INPUT:
% T >> Maturities
% strike >> strike prices (same length of T)
% flag >> enable time dependent parameter upslion

% OUTPUT:
% price 

%--------------------------------------------------------------------------

price = zeros(length(T), length(strike));

% Compute the price for all the maturities and respective strikes

for i = 1:length(T)
    switch(flag)
        case('none')
           params_i = params;
        case('upsilon')
            params_i = [params(1:3), params(i+3)];
        otherwise
            return
    end
    price(i,:) = call_cm_fft(A, n, S0, r(i), T(i), strike, params_i, model);
end

end