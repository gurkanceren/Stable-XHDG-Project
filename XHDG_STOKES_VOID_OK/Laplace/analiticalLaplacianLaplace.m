function lap = analiticalLaplacianLaplace(X)

alpha = 0.1;
beta = 0.3;
a = 5.1;
b = -6.2;
c = 4.3;
d = 3.4;

fun = 'exp(alpha*sin(a*x+b*y)+beta*cos(c*x+d*y))';
dfun_dx = 'eval(fun).*(a*alpha*cos(a*x+b*y)-c*beta*sin(c*x+d*y))';
dfun_dy = 'eval(fun).*(b*alpha*cos(a*x+b*y)-d*beta*sin(c*x+d*y))';
dfun2_dx2 = 'eval(dfun_dx).*(a*alpha*cos(a*x+b*y)-c*beta*sin(c*x+d*y)) - eval(fun).*(a^2*alpha*sin(a*x+b*y)+c^2*beta*cos(c*x+d*y))';
dfun2_dy2 = 'eval(dfun_dy).*(b*alpha*cos(a*x+b*y)-d*beta*sin(c*x+d*y)) - eval(fun).*(b^2*alpha*sin(a*x+b*y)+d^2*beta*cos(c*x+d*y))';

% fun = 'sin(pi*x)*sin(pi*y)';
% dfun2_dx2 = '-pi^2*sin(pi*x)*sin(pi*y)';
% dfun2_dy2 = '-pi^2*sin(pi*x)*sin(pi*y)';

% 
% x = X(:,1);
% y = X(:,2);
% lap = eval(dfun2_dx2) + eval(dfun2_dy2);
% lap = -100*sin(10*x);
% lap = 0;

% dfun2_dx2 = '-pi^2*cos(pi*(x-0.5))*cos(pi*(y-0.5))';
% dfun2_dy2 = '-pi^2*cos(pi*(x-0.5))*cos(pi*(y-0.5))';

x = X(:,1);
y = X(:,2);
lap = eval(dfun2_dx2) + eval(dfun2_dy2);

% lap = zeros(size(lap));