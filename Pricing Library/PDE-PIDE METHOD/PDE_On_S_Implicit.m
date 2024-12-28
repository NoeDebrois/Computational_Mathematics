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

%% PDE NUMERICAL APPROXIMATION:
% ---------- CHANGE HERE IF NOT EU CALL ---------- %
c = max(S' - K, 0); % EU CALL price vector (terminal condition)
% ------------------------------------------------ %

cnew = zeros(size(c)); % new price vector that we are going to fill

% Forward Euler Scheme (Implicit Euler):
A = @(S) - (r * S) * (1 / (2 * dS)) + sigma^2 / 2 * (1 / dS^2) * S^2; % DEPENDS ON S !
B = @(S) (-1 / dt) - sigma^2 * S^2 / dS^2 - r;
C = @(S) (r * S) * (1 / (2 * dS)) + sigma^2 * S^2 / 2 * (1 / dS^2);
% cf page 26 of notes PDE.pdf.

% Build the "big" matrix with A, B, C to solve the systems:
matA = spalloc(N + 1, N + 1, 3 * (N - 1) + 2);
matA(1, 1) = 1;
for row=2:N
    matA(row, [row - 1, row, row + 1]) = [A(S(row)), B(S(row)), C(S(row))];
end
matA(end, end) = 1;

rhs = zeros(N+1, 1);
for j=M:-1:1 % We know c @ t_j -> we compute c @ t_{j-1}...
    % ---------- CHANGE HERE IF NOT EU CALL ---------- %
    rhs(1) = 0; % BC @ Smin
    % ------------------------------------------------ %

    rhs(2:end-1) = - 1/dt * c(2:end-1);

    % ---------- CHANGE HERE IF NOT EU CALL ---------- %
    rhs(end) = Smax - K*exp(-r*(T-(j-1)*dt)); % BC @ Smax = Smin + N*dS
    % ------------------------------------------------ %

    c = matA \ rhs; % Compute c @ t_{j-1} by solving the system
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