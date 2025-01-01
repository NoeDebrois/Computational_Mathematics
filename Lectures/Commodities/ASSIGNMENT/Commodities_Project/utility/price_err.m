function [error] = price_err(vol, F0, K, r, T, A, n, params, model, flag)

% Compute the error between prices computed by Bs and HJM

% OUTPUT
% error >> Difference between mrkt price and model price

%--------------------------------------------------------------------------

mrkt_price = BS(vol, F0, K, r, T);
model_price = HJM(A, n, F0, r, T, K, params, model, flag);

error = model_price - mrkt_price;

end

