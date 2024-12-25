function vel=analyticalVelocityStokes(X,time)

x=X(:,1); y=X(:,2);

%vel = [1+0.*x ; 1+0.*y];
%vel = [x+time ; 0.*y];
%vel = [x+y ; x-y];
%vel = [x.^2 ; -2.*x.*y]; %divergence satisfying vel. field
%vel = [-2.*x.*y ; y.^2]; %divergence satisfying vel. field
%vel=[x.*x.*y ; -x.*y.*y]; %divergence satisfying vel. field

% linear solution:
%vel = [x+y+time ; x-y+time];

% quadratic solution:
vel = [(x.^2)+time ; (-2.*x.*y)+time]; %divergence satisfying vel. field
%vel = [(x.^2)+time ; (-1.*x.*y)+time]; %vel. field not satisfy divergence
%vel = [time.*x ; time.*y]; %vel. field not satisfy divergence
%
%vel = [(-2.*x.*y)+time ; (y.^2)+time]; %divergence satisfying vel. field
%vel = [(-1.*x.*y)+time ; (y.^2)+time]; %vel. field not satisfy divergence
% cubic solution:
%vel=[(x.*x.*y)+time ; (-x.*y.*y)+time]; %divergence satisfying vel. field
%vel=[(x.^3)+time ; (-3*x.*x.*y)+time]; %divergence satisfying vel. field
%vel=[(-3*x.*y.*y)+time ; (y.^3)+time]; %divergence satisfying vel. field
%vel=[(-2*time).*x.*y ; time.*(y.^2)]; %divergence satisfying vel. field

%
%vel = [x+y ; x-y]; %linear solution
%vel = [x.^3 ; -3*(x.^2).*y]; %cubic solution
%vel = [(x.^2).*y ; -(y.^2).*x]; %quadratic solution

%{
Re = 10;
a = Re*Re/4;
b = 4*pi*pi;
lmbda = (Re/2)-sqrt(a+b);
u = 1-exp(2*lmbda.*x).*cos(pi*(4.*y-1));
aa = lmbda/(2*pi);
v = aa*exp(2*lmbda.*x).*sin(pi*(4.*y-1));
vel = [u ; v];
%}

%vel = [y ; x];
%vel = [9+2.*x;15-2.*y];
%vel = [-7.*x+2.*y+5 ; 8.*x+7.*y+11];

%vel = [1.*y ; 0.*x];
