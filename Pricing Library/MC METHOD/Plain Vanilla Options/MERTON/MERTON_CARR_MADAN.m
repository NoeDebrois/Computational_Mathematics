function [Price] = MERTON_CARR_MADAN(K_i, params, T, r, S0)
    Npow = 16;
    N = 2^Npow;
    A = 1000;

    eta = A / N;
    v = [0 : eta : A * (N - 1) / N];
    v(1) = 1e-22;

    lambda = 2 * pi / (N * eta);
    k = -lambda * N / 2 + lambda * (0:N-1);

    CharFunc = @(v) exp(T * CHAR_EXP_MERTON(v, params));
    Z_k = exp(1i * r * v * T) .* (CharFunc(v - 1i) - 1) ./ ...
        (1i * v .* (1i * v + 1));
    
    w = ones(1, N);
    w(1) = 0.5;
    w(end) = 0.5;

    x = w .* eta .* Z_k .* exp(1i * pi * (0:N-1));
    z_k = real(fft(x) / pi);

    C = S0 * (z_k + max(1 - exp(k - r * T), 0));

    K = S0 * exp(k);
    index = find(K > 0.1 * S0 & K < 3 * S0);
    C = C(index);
    K = K(index);

    Price = interp1(K, C, K_i, 'spline');
end
