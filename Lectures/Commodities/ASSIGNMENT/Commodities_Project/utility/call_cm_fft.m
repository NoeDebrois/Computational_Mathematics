function [prices] = call_cm_fft(A, n, S0, r, T, strikes, params, model)

%--------------------------------------------------------------------------

% Compute EU Call price through Carr-Madan formula using FFT for the
% integral numerical computation and kou/merton as jump diffusion model

% INPUT:
% A >> Integral truncation
% n >> 2^n point of integration
% S0 >> Spot price
% r >> rf rate
% T >> maturity
% strikes >> strike prices for the given maturity
% parmas >> model parameters (Kou or Merton)
% model >> 

% OUTPUT:
% prices >> Prices for all the strikes given

%--------------------------------------------------------------------------

% Discretization parameters
% Number of points for the discretization
N = 2^n;
% Discretization step
eta = A/N;
% Discretization points
v = 0:eta:A*(N-1)/N; v(1) = 1e-22;
% Discretization for log-strikes
lambda = 2*pi/(N*eta);
k = -lambda*N/2 + lambda*(0:N-1);

% CharExp of the process Xt
switch(model)
    case('kou')
        [~,~,CharExp] = kou_char_exp(params);
    case('merton')
        [~,~,CharExp] = mertonn_char_exp(params);
    case('VG')
        [~,~,CharExp] = VG_char_exp(params);
    case('NIG')
        [~,~,CharExp] = NIG_char_exp(params);
    case('BS')
        [~,~,CharExp] = BS_char_exp(params);
    otherwise
        print("Invalid model selection");
end

CharFunc = @(u) exp(CharExp(u)*T);
g = exp(1i*r*v*T).*(CharFunc(v - 1i) - 1)./(1i*v.*(1 + 1i*v));

% FFT integral
% Weights (trapezoidal)
w = ones(1,N); w(1) = 0.5; w(end) = 0.5;
x = eta*w.*exp(1i*pi*(0:N-1)).*g;
z_k = real(fft(x)/pi);

% Call prices and strikes
C = S0*(z_k + max(1 - exp(k - r*T),0));
K = S0*exp(k);

% Interpolate on the given strikes
idx = find(K > 0.1*S0 & K < 3*S0);
C = C(idx); K = K(idx);
prices = interp1(K, C, strikes, 'spline');

end

