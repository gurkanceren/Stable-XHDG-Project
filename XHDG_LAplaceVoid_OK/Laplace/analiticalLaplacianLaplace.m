function lap = analiticalLaplacianLaplace(X,c_x,c_y)

%% compute the source term 
%c_x = 1;
%c_y = 1;
 nu = 1;

x = X(:,1);
y = X(:,2);

%{
%***compute du/dx***
trmA = sin(pi*x)+pi.*cos(pi*x);
du_dx = exp(x+y).*sin(pi*y).*trmA;
%***compute d2u/dx2***
trmAA = pi.*cos(pi*x)-pi^2.*sin(pi*x);
trmAB = exp(x+y).*sin(pi*y).*trmAA;
d2u_dx2 = du_dx + trmAB;


%***compute du/dy***
trmB = sin(pi*y)+pi.*cos(pi*y);
du_dy = exp(x+y).*sin(pi*x).*trmB;
%***compute d2u/dy2***
trmBA = pi.*cos(pi*y)-pi^2.*sin(pi*y);
trmBB = exp(x+y).*sin(pi*x).*trmBA;
d2u_dy2 = du_dy + trmBB;
%}


%{
du_dx = 1+0.*x;
du_dy = 0+0.*y;
d2u_dx2 = 0+0.*y;
d2u_dy2 = 0+0.*x;
%}

%{
du_dx = 1+0.*x;
du_dy = 1+0.*y;
d2u_dx2 = 0.*x;
d2u_dy2 = 0.*y;
%}



% du_dx = 0+2.*x;
% du_dy = 0+2.*y;
% d2u_dx2 = 2+0.*x;
% d2u_dy2 = 2+0.*y;



%{
du_dx = 0+1.*y;
du_dy = 0+1.*x;
d2u_dx2 = 0+0.*x;
d2u_dy2 = 0+0.*y;
%}

%{
du_dx = 0+2.*x.*y;
du_dy = 0+x.*x;
d2u_dx2 = 0+2.*y;
d2u_dy2 = 0+0.*y;
%}

%{
du_dx = 0+y.*y;
du_dy = 0+2.*x.*y;
d2u_dx2 = 0+0.*x;
d2u_dy2 = 0+2.*x;
%}

%{
du_dx = 0+3.*x.^2;
du_dy = 0+3.*y.^2;
d2u_dx2 = 0+6.*x;
d2u_dy2 = 0+6.*y;
%}


du_dx = 0+4.*x.^3;
du_dy = 0+4.*y.^3;
d2u_dx2 = 0+12.*x.^2;
d2u_dy2 = 0+12.*y.^2;

%{
du_dx = 0+3.*(x.^2).*y;
du_dy = 0+(x.^3);
d2u_dx2 = 0+6.*x.*y;
d2u_dy2 = 0+0.*y;
%}

%{
du_dx = 0+(y.^3);
du_dy = 0+3.*x.*(y.^2);
d2u_dx2 = 0+0.*x;
d2u_dy2 = 0+6.*x.*y;
%}

lap = -nu.*(d2u_dx2+d2u_dy2); %the source term of the conv-diff eqn.
%disp('Hola')

%lap = lap.*0;



