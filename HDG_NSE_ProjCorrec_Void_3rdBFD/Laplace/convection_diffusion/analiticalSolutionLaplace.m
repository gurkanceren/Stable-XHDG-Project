function u = analiticalSolutionLaplace(X,c_x,c_y)

x = X(:,1);
y = X(:,2);

%c_x = 25;
%c_y = 25;
 nu = 1;

a1 = exp(c_x*(x-1));
a2 = exp(c_y*(y-1));
denm = (1-exp(-c_x))*(1-exp(-c_y));

%% analytical velocity for the steady convection–diffusion case of the Conv-Diff eqn.

u = x.*y.*(1-a1).*(1-a2)/denm;
