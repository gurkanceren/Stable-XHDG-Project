function u_i = vel_initial_condition(X)


x = X(:,1);
y = X(:,2);

%% initial condition for velocity:
%
%% Benchmark example from Ueckermann & Lermusiaux
%u_i = [0.*x ; 0.*y];
%
%
%u_i = [x ; 0.*y];
u_i = [x+y ; x-y]; %linear solution
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


%disp('Hola')



