function vel=analyticalVelocityStokes(X,mu)

x=X(:,1); y=X(:,2);
%vel = [x+y ; x-y]; %linear solution
%vel = [x.^3 ; -3*(x.^2).*y]; %cubic solution
%vel = [(x.^2).*y ; -(y.^2).*x]; %quadratic solution

Re = 10;
a = Re*Re/4;
b = 4*pi*pi;
lmbda = (Re/2)-sqrt(a+b);
u = 1-exp(2*lmbda.*x).*cos(pi*(4.*y-1));
aa = lmbda/(2*pi);
v = aa*exp(2*lmbda.*x).*sin(pi*(4.*y-1));
vel = [u ; v];

%vel = [y ; x];
%vel = [9+2.*x;15-2.*y];
%vel = [-7.*x+2.*y+5 ; 8.*x+7.*y+11];

%vel = [1.*y ; 0.*x];
