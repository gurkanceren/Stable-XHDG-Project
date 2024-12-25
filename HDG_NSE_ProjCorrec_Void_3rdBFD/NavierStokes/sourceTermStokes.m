function s = sourceTermStokes(Xg,Re,time)

x=Xg(:,1); y=Xg(:,2);

%s=Xg*0+2; %linear solution p=x+y+C
%s=Xg*0+0;
%
%s = 0*Xg+0; % constant pressure, p=C
% quadratic solutions:
%s=[0.*x+2-(2/Re) , 0*y+2]; %quadratic solution for p=x+y+C
%
%s=[0*x+0 , 0.*y+0-(2/Re)]; %quadratic solution for p=x+y+C
%
%s=[time+1+0.*x , time+1+0.*y]; %quadratic solution for p=x+y+C
%
% cubic solutions:
%s=[2-(2*y/Re) , 2+(2*x/Re)]; %cubic solution with p=x+y+C
%
%s=[2-(6*x/Re) , 2+(6*y/Re)]; %cubic solution with p=x+y+C
%
%s=[1+time+(6*x/Re) , 1+time-(6*y/Re)]; %cubic solution with p=x+y+C
%
%s=[1-(2.*x.*y) , 1-(2*time/Re)+(y.^2)]; %cubic solution with p=x+y+C

%s=[1-mu*6*x , 1+6*mu*y]; %cubic solution with p=x+y+C
%s=[1-2*mu*y , 1+2*mu*x]; %quadratic solution with p=x+y+C
%
%
%% Benchmark example from Ueckermann & Lermusiaux
%
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
%
%
%{
%% Taylor Vortex problem
a = pi*x;  b = pi*y;  pow = -2*pi*pi*time/Re; 
%
du_dt_x = (2*pi*pi/Re).*cos(a).*sin(b)*exp(pow);
%
dp_dx = (pi/2).*sin(2*a).*exp(2*pow);
%
d2ux_dx2 = (pi*pi).*cos(a).*sin(b).*exp(pow);
%
d2ux_dy2 = (pi*pi).*cos(a).*sin(b).*exp(pow);
%
sx = du_dt_x + dp_dx - (d2ux_dx2+d2ux_dy2)/Re;
%
%
%
du_dt_y = (-2*pi*pi/Re).*sin(a).*cos(b)*exp(pow);
%
dp_dy = (pi/2).*sin(2*b).*exp(2*pow);
%
d2uy_dx2 = (-pi*pi).*sin(a).*cos(b).*exp(pow);
%
d2uy_dy2 = (-pi*pi).*sin(a).*cos(b).*exp(pow);
%
sy = du_dt_y + dp_dy - (d2uy_dx2+d2uy_dy2)/Re;
%
s = [sx , sy];
%}
%
%
%% Kovaszany flow
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