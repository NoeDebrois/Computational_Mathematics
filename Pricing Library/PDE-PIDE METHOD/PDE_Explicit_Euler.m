% Option parameters:
S0 = 1; K = 0.95; r = 0.001; sigma = 0.5; T = 1;

% Space grid:
N = 100; 
Smin = 0.1 * S0;
Smax = 3 * S0;
xmin = log(Smin/S0);
xmax = log(Smax/S0);
dx = (xmax-xmin)/N;
x = xmin + (0:N) * dx; % space grid

% Time grid:
M = 1000;
dt = T/M;
t = (0:M) * dt;

% PDE numerical approximation:
c = max(S0*exp(x') - K, 0); % vector of price
cnew = zeros(size(c)); % new price vector that we are going to fill

A = (r - sigma^2 / 2) * (1 / (2 * dx)) - sigma^2 / 2 * (1 / dx^2);
B = (-1 / dt) + sigma^2 / dx^2 + r;
C = -(r - sigma^2 / 2) * (1 / (2 * dx)) - sigma^2 / 2 * (1 / dx^2);

for j=M:-1:1
    cnew(1) = 0;
    for i=2:N
        cnew(i) = (-dt) * (A * c(i-1) + B * c(i) + C * c(i+1));
    end
    cnew(N+1) = Smax - K*exp(-r*(T-(j-1)*dt));
    c = cnew;
end

figure
plot(S0 * exp(x), c); 
title('Call price at t=0 vs (all possible) Spot Price S_0'); 
xlabel('S at t = 0 (spot price)');
ylabel('Call price at t = 0');

price = interp1(x, c, 0, 'spline') % Interpolation of c at x = 0, i.e S = S0
% (S0 being the one we the one we initialized in "Model parameters").
[price_ex, ~] = blsprice(S0, K, r, T, sigma) % For comparison : B&S price. 