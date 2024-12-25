function p = analyticalpressure(X,time,Re)


x = X(:,1);
y = X(:,2);
C = 0;

%% here constant could be time or a constant.

%% initial condition for pressure:
%p = x+y; %p(x,y,time) = x+y+time, for time=0, p(x,y,time=0)=x+y.
%p = x+y+C;
%
%% Benchmark example from Ueckermann & Lermusiaux
%p = sin(time).*cos(pi.*x).*sin(pi.*y);
%
%
%% Taylor Vortex problem
%{
a = 2*pi*x;  b = 2*pi*y; pow = -4*pi*pi*time/Re;  
%
p = (-0.25).*(cos(a)+cos(b)).*exp(pow); 
%}
%p_i = -x;
%
%% Kovaszany flow
%
%Re = 10;
a = Re*Re/4;
b = 4*pi*pi;
lmbda = (Re/2)-sqrt(a+b);
p = -0.5*exp(2*lmbda*x);
%
%
%disp('Hola')



