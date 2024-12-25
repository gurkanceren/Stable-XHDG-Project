function s = sourceTermStokes(Xg,Re,time)

x=Xg(:,1); y=Xg(:,2);


s=Xg*0+1;
%
%s = 0*Xg+0; % constant pressure, p=C and ux=1, uy=1
%
% linear solutions;
%s=Xg*0+2; %linear solution p=x+y+C
%
% quadratic solutions:
%s=[0.*x+2-(2/Re) , 0*y+2]; %quadratic solution for p=x+y+C
%
%s=[0*x+2 , 0.*y+2-(2/Re)]; %quadratic solution for p=x+y+C
%
%s=[time+1+0.*x , time+1+0.*y]; %quadratic solution for p=x+y+C
%
% cubic solutions:
%s=[2-(2*y/Re) , 2+(2*x/Re)]; %cubic solution with p=x+y+C
%
%s=[2-(6*x/Re) , 2+(6*y/Re)]; %cubic solution with p=x+y+C
%
%s=[2+(6*x/Re) , 2-(6*y/Re)]; %cubic solution with p=x+y+C
%
%s=[1-(2.*x.*y) , 1-(2*time/Re)+(y.^2)]; %cubic solution with p=x+y+C

%s=[1-mu*6*x , 1+6*mu*y]; %cubic solution with p=x+y+C
%s=[1-2*mu*y , 1+2*mu*x]; %quadratic solution with p=x+y+C
%
%
%% Benchmark example from Ueckermann & Lermusiaux
%{
par = 2*pi;
%
cc = 2*(pi^3)/Re;
%
pp = 3*cc;
%
ax = pi.*cos(time).*sin(par.*y).*(sin(pi.*x)).^2;
%
bx = -cc.*sin(time).*sin(par.*y).*(cos(pi.*x)).^2;
%
cx = pp.*sin(time).*sin(par.*y).*(sin(pi.*x)).^2;
%
dx = -pi.*sin(time).*sin(pi.*x).*sin(pi.*y);
%
sx = ax+bx+cx+dx;
%
%
ay = -pi.*cos(time).*sin(par.*x).*(sin(pi.*y)).^2;
%
by = cc.*sin(time).*sin(par.*x).*(cos(pi.*y)).^2;
%
cy = -pp.*sin(time).*sin(par.*x).*(sin(pi.*y)).^2;
%
dy = pi.*sin(time).*cos(pi.*x).*cos(pi.*y);
%
sy = ay+by+cy+dy;
%
s = [sx , sy];
%}
%
%{
Re = 10;
a = Re*Re/4;
b = 4*pi*pi;
lmbda = (Re/2)-sqrt(a+b);
sx = (-2*lmbda)*exp(4*lmbda.*x)+(4*lmbda^2-16*pi^2)*mu*exp(2*lmbda.*x).*cos(pi*(4.*y-1));
sy = (8*pi*lmbda-2*(lmbda^3)/pi)*mu*exp(2*lmbda.*x).*sin(pi*(4.*y-1));
s = [sx , sy];
%}

%s= [1+0*x , 1+0*y]; %%%repeat
%s = [2*x , 0*y];
%s = [3*x.^2, 3*y.^2]; %p=x*x+y*y+C