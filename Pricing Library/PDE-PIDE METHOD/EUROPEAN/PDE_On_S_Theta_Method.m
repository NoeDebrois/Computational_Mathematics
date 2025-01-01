%% PARAMETERS:
% Option parameters:
S0 = 1; K = 0.95; r = 0.001; sigma = 0.5; T = 1;

% Space grid:
N = 1000; 
Smin = 0.1 * S0;
Smax = 3 * S0;
dS = (Smax-Smin)/N;
S = Smin + (0:N) * dS; % space grid (NO LOG PRICE! It's the real price S)

% Time grid:
M = 1000;
dt = T/M;
t = (0:M) * dt;

% Theta method:
theta = 0.5;
% 0 : Implicit Euler, 0.5 : Crank-Nicholson, 1 : Explicit Euler

%% PDE NUMERICAL APPROXIMATION:
% ---------- CHANGE HERE IF NOT EU CALL ---------- %
c = max(S' - K, 0); % EU CALL price vector (terminal condition)
% ------------------------------------------------ %
cnew = zeros(size(c)); % new price vector that we are going to fill

% Build M1 (cf my notes):
A = @(S) (1 - theta) * dt * (-r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2)); % DEPENDS ON S !
B = @(S) -1 + dt * (1 - theta) * (-sigma^2 * S^2 / (dS^2) - r);
C = @(S) dt * (1 - theta) * (r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));
% cf page 26 of notes PDE.pdf.

Mat = spalloc(N + 1, N + 1, 3 * (N - 1) + 2); % M1 in my notes
Mat(1, 1) = 1; 
for row=2:N
    Mat(row, [row - 1, row, row + 1]) = [A(S(row)) B(S(row)) C(S(row))];
end
Mat(end, end) = 1;

% Build M2 (cf my notes):
Ah = @(S) -(theta) * dt * (-r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2)); % DEPENDS ON S !
Bh = @(S) -1 - dt * (theta) * (-sigma^2 * S^2 / (dS^2) - r);
Ch = @(S) -dt * (theta) * (r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));

Mat_rhs = spalloc(N + 1, N + 1, 3 * (N - 1)); % M2 in my notes
for row=2:N
    Mat_rhs(row, [row - 1, row, row + 1]) = [Ah(S(row)) Bh(S(row)) Ch(S(row))];
end

rhs = zeros(N+1, 1);
for j=M:-1:1 % We know c @ t_j -> we compute c @ t_{j-1}...
    rhs = Mat_rhs * c; % M2 * V_{j+1}

    % ---------- CHANGE HERE IF NOT EU CALL ---------- %
    rhs(1) = 0; % BC @ Smin
    % ------------------------------------------------ %

    % ---------- CHANGE HERE IF NOT EU CALL ---------- %
    rhs(end) = Smax - K*exp(-r*(T-(j-1)*dt)); % BC @ Smax = Smin + N*dS
    % ------------------------------------------------ %
    
    c = Mat \ rhs; % M1^{-1} * (M2 * V_{j+1} + BC_j)
end

%% PLOT:
figure
plot(S, c); 
title('Call price at t=0 vs (all possible) Spot Price S_0'); 
xlabel('S at t = 0 (spot price)');
ylabel('Call price at t = 0');

%% INTERPOLATION:
price_PDE = interp1(S, c, S0, 'spline') % Interpolation of c at spot price S0 (fixed above).
[price_CHECK, ~] = blsprice(S0, K, r, T, sigma) % Comparison w/ B&S price. 