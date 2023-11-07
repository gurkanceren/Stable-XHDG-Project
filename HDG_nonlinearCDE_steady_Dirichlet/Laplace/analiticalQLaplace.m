function q = analiticalQLaplace(X)
global mu
%% compute the source term for the conv-diff eqn.

x = X(:,1);
y = X(:,2);

q=zeros(2*size(x,1),1);

% a1 = exp(c_x*(x-1));
% a2 = exp(c_y*(y-1));
% a = 1-exp(-c_x);
% b = 1-exp(-c_y);
% denm = a*b;

%% analytical velocity for the steady convection–diffusion case of the Conv-Diff eqn.
%u = x*y*(1-a1)*(1-a2)/denm;

%***compute du/dx***
% bb = 1+(x*c_x);
% b_n = y.*(1-bb.*a1).*(1-a2);
% du_dx = b_n/denm;
du_dx = 0+x.*2; 
q_x = -mu*du_dx;

%***compute du/dy***
% aa = 1+(y*c_y);
% a_n = x.*(1-a1).*(1-aa.*a2);
% du_dy = a_n/denm;
du_dy = 0+y.*0; 
q_y = -mu*du_dy;

%% ***compute qt***
%qt_x = q_x + c_x*u;
%qt_y = q_y + c_y*u;

q(1:2:end) = q_x; 
q(2:2:end) = q_y; 




