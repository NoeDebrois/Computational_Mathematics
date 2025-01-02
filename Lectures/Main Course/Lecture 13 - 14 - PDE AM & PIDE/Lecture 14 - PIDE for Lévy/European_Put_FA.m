clear
close all
% Price a European Put option under the Merton model
% PIDE -> logprice, finite activity, implicit scheme (operator splitting)

%% Input & grids
S0=1; K=1; T=0.5; r=0.1/100; sigma=0.4; lambda=2; muJ=-0.01; deltaJ=0.2;
M=250; dt=T/M;
N=1000;
Smin=0.3*S0; Smax=3*S0; 
xmin=log(Smin/S0); xmax=log(Smax/S0); 
dx=(xmax-xmin)/N; x=xmin+(0:N)'*dx;

%% Deal with the integral (truncation, lambda, alpha)
k=@(y) lambda*exp(-(y-muJ).^2/(2*deltaJ^2))/(deltaJ*sqrt(2*pi));
tol=1e-8;
tmin=-0.5;
while (abs(k(tmin))>tol)
    tmin=tmin-0.5;
end
tmax=0.5;
while (abs(k(tmax))>tol)
    tmax=tmax+0.5;
end
figure
tt=linspace(tmin,tmax,2*N);
plot(tt,k(tt));
lambda_num=trapz(tt,k(tt))
alpha=trapz(tt,(exp(tt)-1).*k(tt))

%% Deal with the implicit matrix
Mat=spalloc(N+1,N+1,3*(N-1)+2);
Mat(1,1)=1; Mat(end,end)=1;
A=-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2);
B=-1/dt-sigma^2/(dx^2)-(r+lambda);
C= (r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2);
for i=2:N
    Mat(i,[i-1 i i+1])=[A B C];
end

%% Backward in time
c=max(0,K-S0*exp(x));
rhs=zeros(N+1,1); 
for j=M:-1:1
    rhs(1)=K*exp(-r*(T-(j-1)*dt))-S0*exp(xmin); % BC
    rhs(end)=0; %BC
    I=integral_levy(x,c,tt,k,K*exp(-r*(T-(j-1)*dt)),S0);
    rhs(2:end-1)=-1/dt*c(2:end-1)-I;
    c=Mat\rhs;
end
figure
plot(S0*exp(x),c);
price=interp1(x,c,0,'spline')
price_cm_eurcall=FFT_CM_Call_Merton(K,[sigma lambda muJ deltaJ],T,r,S0);
price_cm_eurput=price_cm_eurcall-S0+K*exp(-r*T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra function
function I=integral_levy(x,c,tt,k,Kdisc,S0)

I=zeros(length(x)-2,1);
for i=2:length(x)-1
    I(i-1)=trapz(tt, cfun(x(i)+tt,x,c,Kdisc,S0).*k(tt) );
end
end

function c_int=cfun(y,x,c,Kdisc,S0)

c_int=zeros(size(y));
index=find(y<=x(1));
c_int(index)=Kdisc-S0*exp(y(index));
%index=find(y>=x(end));
%c_int(index)=0;
index=find( (y>x(1)).*(y<x(end)) );
c_int(index)=interp1(x,c,y(index));

end