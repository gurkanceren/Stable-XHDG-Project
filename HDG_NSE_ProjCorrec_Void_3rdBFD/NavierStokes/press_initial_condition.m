function p_i = press_initial_condition(X,time)


x = X(:,1);
y = X(:,2);
C = 0;

%% here constant could be time or a constant.

%% initial condition for pressure:
%p_i = x+y; %p(x,y,time) = x+y+time, for time=0, p(x,y,time=0)=x+y.
%
%% Benchmark example from Ueckermann & Lermusiaux
p_i = C.*(x+y);
%
%{
%% Taylor Vortex problem
a = 2*pi*x;  b = 2*pi*y;  
%
p_i = (-0.25).*(cos(a)+cos(b)); 
%}
%
%
%p_i = -x;

%disp('Hola')



