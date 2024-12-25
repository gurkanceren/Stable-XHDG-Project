function u0 = analiticalSolutionLaplace(X)

alpha = 0.1;
beta = 0.3;
a = 5.1;
b = -6.2;
c = 4.3;
d = 3.4;

x = X(:,1);
y = X(:,2);

%u0=x;
u0 = exp(alpha*sin(a*x+b*y)+beta*cos(c*x+d*y));

%u0=x.^2;

 %u0 = sin(pi*x).*sin(pi*y);
 
% u0 = exp(x);

% u0 = cos(pi*(x-0.5)).*cos(pi*(y-0.5));


% u0 = a*x + b*y;