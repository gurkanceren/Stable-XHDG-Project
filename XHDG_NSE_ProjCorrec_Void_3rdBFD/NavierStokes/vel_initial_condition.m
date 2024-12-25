function u_i = vel_initial_condition(X)

x = X(:,1);
y = X(:,2);

%% initial condition for velocity:
%
%% Benchmark example from Ueckermann & Lermusiaux
%u_i = [0.*x ; 0.*y];
%
%
%% Taylor Vortex problem
%{
a = pi*x;  b = pi*y;  
%
u_x = -1.*cos(a).*sin(b);
%
u_y = sin(a).*cos(b);
%
u_i = [u_x ; u_y];
%}
%
%
%u_i = [x+y ; x-y]; %linear solution
% quadratic solution:
%u_i = [x.^2 ; -2.*x.*y];
%u_i = [x.^2 ; -1.*x.*y];
%u_i = [1+0.*x ; 1+0.*y];
%
%u_i = [-2.*x.*y ; y.^2];
%u_i = [-1.*x.*y ; y.^2];
% cubic solution:
%u_i = [(x.^2).*y ; -(y.^2).*x];
%u_i = [x.^3 ; -3*(x.^2).*y];
%u_i = [-3*x.*y.*y ; y.^3];
%u_i = [0.*x ; 0.*y];

%u_i = [1+0.*x ; 1+0.*y];
%
%% Kovaszany flow
%
Re = 10;
a = Re*Re/4;
b = 4*pi*pi;
lmbda = (Re/2)-sqrt(a+b);
u = 1-exp(2*lmbda.*x).*cos(pi*(4.*y-1));
aa = lmbda/(2*pi);
v = aa*exp(2*lmbda.*x).*sin(pi*(4.*y-1));
vel = [u ; v];
%
%

%disp('Hola')



