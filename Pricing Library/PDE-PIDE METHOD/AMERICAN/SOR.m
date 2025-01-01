function xnew = SOR(A, rhs, x)
% Initialization:
xnew = zeros(size(x)); 
n = length(A);

% Solver parameters:
maxiter = 100;
tol = 1e-7; 
w = 1.5; % "learning rate"

% PSOR:
for i=1:maxiter
    for j=1:n % Size of the vector
        % Implementation of the formula (cf page 39):
        y = (rhs(j) - A(j, 1:j-1) * xnew(1:j-1, 1) ...
            - A(j, j+1:end) * x(j+1:end, 1)) / A(j, j);
        xnew(j) = x(j) + w * (y - x(j));
    end
    if norm(xnew - x) < tol
        break
    else
        x = xnew;
    end
end