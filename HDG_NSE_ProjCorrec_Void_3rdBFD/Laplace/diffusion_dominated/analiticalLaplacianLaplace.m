function lap = analiticalLaplacianLaplace(X,c_x,c_y)

%{
alpha = 0.1;
beta = 0.3;
a = 5.1;
b = -6.2;
c = 4.3;
d = 3.4;

x = X(:,1);
y = X(:,2);

lap=x*0;
A=0.2;
 
% fun = 'exp(alpha*sin(a*x+b*y)+beta*cos(c*x+d*y))';
% dfun_dx = 'eval(fun).*(a*alpha*cos(a*x+b*y)-c*beta*sin(c*x+d*y))';
% dfun_dy = 'eval(fun).*(b*alpha*cos(a*x+b*y)-d*beta*sin(c*x+d*y))';
% dfun2_dx2 = 'eval(dfun_dx).*(a*alpha*cos(a*x+b*y)-c*beta*sin(c*x+d*y)) - eval(fun).*(a^2*alpha*sin(a*x+b*y)+c^2*beta*cos(c*x+d*y))';
% dfun2_dy2 = 'eval(dfun_dy).*(b*alpha*cos(a*x+b*y)-d*beta*sin(c*x+d*y)) - eval(fun).*(b^2*alpha*sin(a*x+b*y)+d^2*beta*cos(c*x+d*y))';

% fun = 'sin(pi*x)*sin(pi*y)';
% dfun2_dx2 = '-pi^2.*sin(pi*x).*sin(pi*y)';
% dfun2_dy2 = '-pi^2.*sin(pi*x).*sin(pi*y)';

% 
% x = X(:,1);
% y = X(:,2);
% lap = eval(dfun2_dx2) + eval(dfun2_dy2);
% lap = -100*sin(10*x);
% lap = 0;

% dfun2_dx2 = '-pi^2*cos(pi*(x-0.5))*cos(pi*(y-0.5))';
% dfun2_dy2 = '-pi^2*cos(pi*(x-0.5))*cos(pi*(y-0.5))';

for i=1:length(lap)
    
    %r(i)=sqrt(x(i)^2+y(i)^2);
    
    if   x(i)<0.2  % r(i)>0.4     % Domain 1
     
        
    lap=100*x^3;    
    %lap= -2*exp(sin(x^2))*[(2*x^2*sin(x^2)-2*x^2*(cos(x^2))^2-cos(x^2))];
        

    else
      
    lap=100*x^3;    
    %lap= -2*exp(sin(x^2))*[(2*x^2*sin(x^2)-2*x^2*(cos(x^2))^2-cos(x^2))];          
    %lap= A*-2*exp(sin(x^2))*[(2*x^2*sin(x^2)-2*x^2*(cos(x^2))^2-cos(x^2))];     
  
    end
    
end
%lap = eval(dfun2_dx2) + eval(dfun2_dy2);

%lap = zeros(size(lap));

%lap=-pi^2.*sin(pi*x);

%lap = -4*pi^2.*sin(pi*x).*cos(pi*x);
%}

%% compute the source term for the conv-diff eqn. (diffusion dominated problem).
%c_x = 1;
%c_y = 1;
 nu = 1;

x = X(:,1);
y = X(:,2);

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

%{
du_dx = 0+2.*x;
du_dy = 0+2.*y;
d2u_dx2 = 2+0.*x;
d2u_dy2 = 2+0.*y;
%}

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

%{
du_dx = 0+4.*x.^3;
du_dy = 0+4.*y.^3;
d2u_dx2 = 0+12.*x.^2;
d2u_dy2 = 0+12.*y.^2;
%}

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

lap = c_x.*du_dx+c_y.*du_dy-nu.*(d2u_dx2+d2u_dy2); %the source term of the conv-diff eqn.
%disp('Hola')

%lap = 1;



