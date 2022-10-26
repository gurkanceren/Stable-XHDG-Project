function qx = grad_ux(X)

alpha = 0.1;
beta = 0.3;
a = 5.1;
b = -6.2;
c = 4.3;
d = 3.4;

fun = 'exp(alpha*sin(a*x+b*y)+beta*cos(c*x+d*y))';
dfun_dx = 'eval(fun).*(a*alpha*cos(a*x+b*y)-c*beta*sin(c*x+d*y))';

x = X(:,1);
y = X(:,2);
qx = eval(dfun_dx) ;