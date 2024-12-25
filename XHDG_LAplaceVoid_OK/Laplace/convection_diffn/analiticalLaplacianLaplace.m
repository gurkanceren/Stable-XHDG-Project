function lap = analiticalLaplacianLaplace(X,c_x,c_y)

%% compute the source term for the conv-diff eqn. (diffusion dominated problem).
%c_x = 25;
%c_y = 25;
 nu = 1;

x = X(:,1);
y = X(:,2);

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

%***compute d2u/dx2***
trmAB = -c_x.*y.*a1.*(x*c_x+2).*(1-a2); 
d2u_dx2 = trmAB/denm;

%***compute d2u/dy2***
trmBA = -c_y.*x.*(y*c_y+2).*a2.*(1-a1);
d2u_dy2 = trmBA/denm;

%*****the source term of the conv-diff eqn.*****
lap = c_x.*du_dx+c_y.*du_dy-nu.*(d2u_dx2+d2u_dy2); 

%disp('Hola')

