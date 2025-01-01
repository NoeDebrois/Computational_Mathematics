clear
close all
% Pricing a plain vanilla Put Option
% B&S model - logprice PDE
% ThetaMethod
%% 1. Parameters & Grids
S0=1; K=0.95; r=0.001; sigma=0.5; T=1;
M=100;
N=500;
theta=0.5; % 0 Implicit, 0.5 CN, 1 Explicit
Smin=0.1*S0; Smax=3*S0;
xmin=log(Smin/S0); xmax=log(Smax/S0);
dx=(xmax-xmin)/N;
%x=linspace(xmin,xmax,N+1);
x=xmin+(0:N)*dx;
dt=T/M;
t=(0:M)*dt;

%% Matrix Construction
A=(1-theta)*dt*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
B=-1+dt*(1-theta)*(-sigma^2/(dx^2)-r);
C=dt*(1-theta)*( (r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
Mat=spalloc(N+1,N+1,3*(N-1)+2);
Mat(1,1)=1;
for i=2:N
    Mat(i, [i-1,i,i+1])=[A B C];
end
Mat(end,end)=1;
Ah=-(theta)*dt*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
Bh=-1-dt*(theta)*(-sigma^2/(dx^2)-r);
Ch=-dt*(theta)*( (r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
Mat_rhs=spalloc(N+1,N+1,3*(N-1));
for i=2:N
    Mat_rhs(i, [i-1,i,i+1])=[Ah Bh Ch];
end

%% Backward in time loop
p=max(K-S0*exp(x'),0); % Payoff at maturity
rhs=zeros(N+1,1);
for j=M:-1:1
    % know c t_j -> compute cnew t_{j-1}
    rhs=Mat_rhs*p; % M2 * [V_{j+1},i] cf page 22
    rhs(1)=K*exp(-r*(T-(j-1)*dt))-Smin; % BC_j (cf page 22 except that here it's a PUT)
    rhs(end)=0; % PUT 
    p=SOR(Mat,rhs,p); % Slower than p=Mat\rhs;
    % vs notes page 22 : Mat=M1, rhs=M2 * [V_{j+1},i] + BCj.
end
figure
plot(S0*exp(x),p); title('put price'); xlabel('S');
price=interp1(x,p,0,'spline')





[~,price_ex]=blsprice(S0,K,r,T,sigma)