function u = analiticalSolutionLaplace(X,time)

x = X(:,1);
y = X(:,2);

%u = 0+x.*time;
%u = (x+y).*time;
%u = (x.^2+y.^2).*time;
%u = x.*y.*time;
%u = (x.^2).*y.*time;
%u = x.*(y.^2).*time;
%u = (x.^3+y.^3).*time;
%u = (x.^4+y.^4).*time;
%u = (x.^3).*y.*time;
%u = x.*(y.^3).*time;

%************************************************************************

%u = 0+1.*x;
%u = x+y;
%u = x.^2+y.^2;
%u = (x.*y);
%u = ((x.^2).*y);
%u = x.*(y.^2);
%u = (x.^3+y.^3);
%u = (x.^4+y.^4);
%u = ((x.^3).*y);
u = x.*(y.^3);
%u = (x.^5+y.^5)+10*time;
%u = (x.^6+y.^6)+10*time;
%{
u_0 = 0.8; v_0 = 0.8;
alpha_x = 0.01; alpha_y = 0.01;
x_0 = 0.5; y_0 = 0.5;

a = 1/(1+4*time);

bb_1 = x-(u_0*time)-0.5;
bb_2 = (1+4*time)*alpha_x;
bb_3 = (bb_1).^2;
b = -bb_3/bb_2;

cc_1 = y-(v_0*time)-0.5;
cc_2 = (1+4*time)*alpha_y;
cc_3 = (cc_1).^2;
c = -cc_3/cc_2;

u = a.*exp(b).*exp(c);
%}

%disp('Hola');
