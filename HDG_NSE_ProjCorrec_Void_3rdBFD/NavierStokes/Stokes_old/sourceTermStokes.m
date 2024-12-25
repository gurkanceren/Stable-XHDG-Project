function s = sourceTermStokes(Xg,mu)

x=Xg(:,1); y=Xg(:,2);

%s=Xg*0+1; %linear solution p=x+y+C
%s=[1-mu*6*x , 1+6*mu*y]; %cubic solution with p=x+y+C
%s=[1-2*mu*y , 1+2*mu*x]; %quadratic solution with p=x+y+C

Re = 10;
a = Re*Re/4;
b = 4*pi*pi;
lmbda = (Re/2)-sqrt(a+b);
sx = (-2*lmbda)*exp(4*lmbda.*x)+(4*lmbda^2-16*pi^2)*mu*exp(2*lmbda.*x).*cos(pi*(4.*y-1));
sy = (8*pi*lmbda-2*(lmbda^3)/pi)*mu*exp(2*lmbda.*x).*sin(pi*(4.*y-1));
s = [sx , sy];

%s= [1+0*x , 1+0*y]; %%%repeat
%s = [2*x , 0*y];
%s = [3*x.^2, 3*y.^2]; %p=x*x+y*y+C