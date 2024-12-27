%% PARAMETERS:
% Option parameters:
S0 = 1; K = 0.95; r = 0.001; sigma = 0.5; T = 1;

% Space grid:
N = 1000; 
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

% Theta method:
theta = 0.5;
% 0 : Implicit Euler, 0.5 : Crank-Nicholson, 1 : Explicit Euler

%% PDE NUMERICAL APPROXIMATION:
% ---------- CHANGE HERE IF NOT EU CALL ---------- %
c = max(S0*exp(x') - K, 0); % EU CALL price vector (terminal condition)
% ------------------------------------------------ %
cnew = zeros(size(c)); % new price vector that we are going to fill

% Build M1 (cf my notes):
A = (1 - theta) * dt * (-(r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));
B = -1 + dt * (1 - theta) * (-sigma^2 / (dx^2) - r);
C = dt * (1 - theta) * ((r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));

Mat = spalloc(N + 1, N + 1, 3 * (N - 1) + 2); % M1 in my notes
Mat(1, 1) = 1; 
for i=2:N
    Mat(i, [i - 1, i, i + 1]) = [A B C];
end
Mat(end, end) = 1;

% Build M2 (cf my notes):
Ah = - (theta) * dt * (-(r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));
Bh = -1 - dt * (theta) * (-sigma^2 / (dx^2) - r);
Ch = - dt * (theta) * ((r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));

Mat_rhs = spalloc(N + 1, N + 1, 3 * (N - 1)); % M2 in my notes
for i=2:N
    Mat_rhs(i, [i - 1, i, i + 1]) = [Ah Bh Ch];
end

rhs = zeros(N+1, 1);
for j=M:-1:1 % We know c @ t_j -> we compute c @ t_{j-1}...
    rhs = Mat_rhs * c; % M2 * V_{j+1}

    % ---------- CHANGE HERE IF NOT EU CALL ---------- %
    rhs(1) = 0; % BC @ xmin = 0
    % ------------------------------------------------ %

    % ---------- CHANGE HERE IF NOT EU CALL ---------- %
    rhs(end) = Smax - K*exp(-r*(T-(j-1)*dt)); % BC @ xmax = xmin + N*dx
    % ------------------------------------------------ %
    
    c = Mat \ rhs; % M1^{-1} * (M2 * V_{j+1} + BC_j)
end

%% PLOT:
figure
plot(S0 * exp(x), c); 
title('Call price at t=0 vs (all possible) Spot Price S_0'); 
xlabel('S at t = 0 (spot price)');
ylabel('Call price at t = 0');

%% INTERPOLATION:
price_PDE = interp1(x, c, 0, 'spline') % Interpolation of c at x = 0, i.e S = S0
% (S0 being the one we the one we initialized in "Model parameters").
[price_CHECK, ~] = blsprice(S0, K, r, T, sigma) % Comparison w/ B&S price. 