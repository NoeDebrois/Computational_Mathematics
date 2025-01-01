% A basic implementation of SOR. But when we'll use it, it will be a
% function (cf SOR.m).

clear
A=[ 2 1 0.2; 0.1 3 0.2; 3 5 12];
rhs=[1;3;4];
x=A\rhs

% SOR
x=zeros(3,1); % guess
xnew=zeros(3,1); 
n=length(A);
maxiter=100;
tol=1e-7;
w=1.5;
for i=1:maxiter
    for j=1:n
        y=(rhs(j)-A(j,1:j-1)*xnew(1:j-1,1)...
            -A(j,j+1:end)*x(j+1:end,1))/A(j,j);
        xnew(j)=x(j)+w*(y-x(j));
    end
    if norm(xnew-x)<tol
        break
    else
        x=xnew;
    end
end
xnew   
xnewSOR = SOR(A, rhs, x)