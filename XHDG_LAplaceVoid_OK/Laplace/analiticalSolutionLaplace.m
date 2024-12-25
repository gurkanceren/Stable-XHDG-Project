function u = analiticalSolutionLaplace(X)


x = X(:,1);
y = X(:,2);


%% analytical velocity for the diffusion dominated case of the Conv-Diff eqn.

%u = exp(x+y).*sin(pi*x).*sin(pi*y);

%u = 0+1.*x; 
%u = x+y;
%u = x.^2+y.^2;
%u = x.*y;
%u = (x.^2).*y;
%u = x.*(y.^2);
%u = x.^3+y.^3;
u = x.^4+y.^4;
%u = (x.^3).*y;
%u = x.*(y.^3);