function I = INTEGRAL_LEVY(x, V, ynodes, nu, lb_BC, ub_BC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integral part of the numerical approximation of option pricess via PIDE 
% on logprices in Lévy framework for General Levy processes
% INPUT:
% - x  = logprices grid
% - V  = last evaluation of the option price
% - ynodes = integration nodes
% - nu = function, Lévy measure 
% - lb_BC & ub_BC = BC for values y in ynodes outside of [xmin, xmax]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = zeros(length(x) - 2, 1);
dx = x(2) - x(1);

for i=2:length(x)-1
    I(i - 1) = trapz(ynodes, (Vfun(x(i) + ynodes, x, V, lb_BC, ub_BC) ...
               - V(i) - (exp(ynodes) - 1) * (V(i + 1) - V(i - 1)) / (2 * dx)) .* nu(ynodes));
end
end

function V_int = Vfun(y, x, V, lb_BC, ub_BC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper: computes the value of V on the grid y adjusting BC.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_int = zeros(size(y));

% Values of V for y in [xmin, xmax]:
index = find((y > x(1)) .* (y < x(end))); % y INSIDE [xmin, xmax] !
V_int(index) = interp1(x, V, y(index));   % Interpolation of V.

% Values of V for y < xmin:
index = find(y < x(1));
V_int(index) = lb_BC(y(index));

% Values of V for y > xmax:
index = find(y > x(end));
V_int(index) = ub_BC(y(index));
end