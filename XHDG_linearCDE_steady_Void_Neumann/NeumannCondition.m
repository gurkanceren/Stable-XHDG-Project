
function [ge]=NeumannCondition(n_g,X,c_x,c_y,mu)
% delu/deln=gradu*n
x = X(:,1); y = X(:,2);

%du_dx = 1+0.*x;
%du_dy = 1+0.*y;

u = analiticalSolutionLaplace(X);

%
%***compute du/dx***
trmA = sin(pi*x)+pi.*cos(pi*x);
du_dx = exp(x+y).*sin(pi*y).*trmA;
%***compute du/dy***
trmB = sin(pi*y)+pi.*cos(pi*y);
du_dy = exp(x+y).*sin(pi*x).*trmB;
%*******************
%

%{
a1 = exp(c_x*(x-1));
a2 = exp(c_y*(y-1));
denm = (1-exp(-c_x))*(1-exp(-c_y));

%***compute du/dx***
bb = 1+(x*c_x);
b_n = y.*(1-bb.*a1).*(1-a2);
du_dx = b_n/denm;

%***compute du/dy***
aa = 1+(y*c_y);
a_n = x.*(1-a1).*(1-aa.*a2);
du_dy = a_n/denm;
%*******************
%}

gradu=[du_dx;du_dy];
q = -mu.*gradu;

c = [c_x;c_y];
cu = c*u;

aa = cu+q;
ge=n_g*aa;

%disp('Hola');




