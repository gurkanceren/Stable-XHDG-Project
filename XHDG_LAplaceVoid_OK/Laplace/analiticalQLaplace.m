function q = analiticalQLaplace(X,c_x,c_y)


%% compute the source term for the conv-diff eqn.
c_x = 25;
c_y = 25;
 nu = 1;
visc = 1;

x = X(:,1);
y = X(:,2);

a1 = exp(c_x*(x-1));
a2 = exp(c_y*(y-1));
a = 1-exp(-c_x);
b = 1-exp(-c_y);
denm = a*b;

%% analytical velocity for the steady convection–diffusion case of the Conv-Diff eqn.
u = x*y*(1-a1)*(1-a2)/denm;

%***compute du/dx***
bb = 1+(x*c_x);
b_n = y.*(1-bb.*a1).*(1-a2);
du_dx = b_n/denm;
q_x = -visc*du_dx;

%***compute du/dy***
aa = 1+(y*c_y);
a_n = x.*(1-a1).*(1-aa.*a2);
du_dy = a_n/denm;
q_y = -visc*du_dy;

%% ***compute qt***
%qt_x = q_x + c_x*u;
%qt_y = q_y + c_y*u;

q = [q_x ; q_y];




