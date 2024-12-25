function p_i = press_initial_condition(X,time)


x = X(:,1);
y = X(:,2);
C = 0;

%% here constant could be time or a constant.

%% initial condition for pressure:
%p_i = x+y; %p(x,y,time) = x+y+time, for time=0, p(x,y,time=0)=x+y.
%
%p_i = C.*(x+y);
%
%% Benchmark example from Ueckermann & Lermusiaux
%
%p_i = sin(C).*cos(pi.*x).*sin(pi.*y);
%
%
%% Taylor Vortex problem
%{
a = 2*pi*x;  b = 2*pi*y;  
%
p_i = (-0.25).*(cos(a)+cos(b)); 
%}
%
%% Kovaszany flow
%
Re = 10;
a = Re*Re/4;
b = 4*pi*pi;
lmbda = (Re/2)-sqrt(a+b);
p = -0.5*exp(2*lmbda*x);
%
%
%p_i = -x;

%disp('Hola')



