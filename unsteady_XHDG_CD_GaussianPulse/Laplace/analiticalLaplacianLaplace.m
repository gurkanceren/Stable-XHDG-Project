function lap = analiticalLaplacianLaplace(X,c_x,c_y,time)


%% compute the source term for the conv-diff eqn. (diffusion dominated problem).
%c_x = 1;
%c_y = 1;
 nu = 0.01;

x = X(:,1);
y = X(:,2);

%{
du_dt = x;
du_dx = time+0.*x;
du_dy = 0+0.*y;
d2u_dx2 = 0+0.*y;
d2u_dy2 = 0+0.*x;
%}

%{
du_dt = x+y;
du_dx = time+0.*x;
du_dy = time+0.*y;
d2u_dx2 = 0+0.*y;
d2u_dy2 = 0+0.*x;
%}


%{
du_dt = x.^2+y.^2;
du_dx = 0+2.*x.*time;
du_dy = 0+2.*y.*time;
d2u_dx2 = 2*time+0.*x;
d2u_dy2 = 2*time+0.*y;
%}

%{
du_dt = x.*y;
du_dx = 0+1.*y*time;
du_dy = 0+1.*x*time;
d2u_dx2 = 0+0.*x;
d2u_dy2 = 0+0.*y;
%}

%{
du_dt = (x.^2).*y;
du_dx = 0+2.*x.*y.*time;
du_dy = 0+x.*x.*time;
d2u_dx2 = 0+2.*y.*time;
d2u_dy2 = 0+0.*y;
%}

%{
du_dt = x.*(y.^2);
du_dx = 0+y.*y.*time;
du_dy = 0+2.*x.*y.*time;
d2u_dx2 = 0+0.*x;
d2u_dy2 = 0+2.*x.*time;
%}

%{
du_dt = x.^3+y.^3;
du_dx = 0+3.*x.^2.*time;
du_dy = 0+3.*y.^2.*time;
d2u_dx2 = 0+6.*x.*time;
d2u_dy2 = 0+6.*y.*time;
%}

%{
du_dt = x.^4+y.^4;
du_dx = 0+4.*x.^3.*time;
du_dy = 0+4.*y.^3.*time;
d2u_dx2 = 0+12.*x.^2.*time;
d2u_dy2 = 0+12.*y.^2.*time;
%}

%{
du_dt = (x.^3).*y;
du_dx = 0+3.*(x.^2).*y.*time;
du_dy = 0+(x.^3).*time;
d2u_dx2 = 0+6.*x.*y.*time;
d2u_dy2 = 0+0.*y;
%}

%{
du_dt = x.*(y.^3);
du_dx = 0+(y.^3).*time;
du_dy = 0+3.*x.*(y.^2).*time;
d2u_dx2 = 0+0.*x;
d2u_dy2 = 0+6.*x.*y.*time;
%}

%**************************************************************

%{
du_dt = 1+0.*x;
du_dx = 1+0.*x;
du_dy = 0+0.*y;
d2u_dx2 = 0+0.*y;
d2u_dy2 = 0+0.*x;
%}

%{
du_dt = 1+0.*x;
du_dx = 1+0.*x;
du_dy = 1+0.*y;
d2u_dx2 = 0+0.*y;
d2u_dy2 = 0+0.*x;
%}

%{
du_dt = 1+0.*x;
du_dx = 0+2.*x;
du_dy = 0+2.*y;
d2u_dx2 = 2+0.*x;
d2u_dy2 = 2+0.*y;
%}

%{
du_dt = 1+0.*x;
du_dx = 0+1.*y;
du_dy = 0+1.*x;
d2u_dx2 = 0+0.*x;
d2u_dy2 = 0+0.*y;
%}

%{
du_dt = 1+0.*x;
du_dx = 0+2.*x.*y;
du_dy = 0+x.*x;
d2u_dx2 = 0+2.*y;
d2u_dy2 = 0+0.*y;
%}

%{
du_dt = 1+0.*x;
du_dx = 0+y.*y;
du_dy = 0+2.*x.*y;
d2u_dx2 = 0+0.*x;
d2u_dy2 = 0+2.*x;
%}

%{
du_dt = 1+0.*x;
du_dx = 0+3.*x.^2;
du_dy = 0+3.*y.^2;
d2u_dx2 = 0+6.*x;
d2u_dy2 = 0+6.*y;
%}

%{
du_dt = 10+0.*x;
du_dx = 0+4.*x.^3;
du_dy = 0+4.*y.^3;
d2u_dx2 = 0+12.*x.^2;
d2u_dy2 = 0+12.*y.^2;
%}

%{
du_dt = 1+0.*x;
du_dx = 0+3.*(x.^2).*y;
du_dy = 0+(x.^3);
d2u_dx2 = 0+6.*x.*y;
d2u_dy2 = 0+0.*y;
%}

%{
du_dt = 1+0.*x;
du_dx = 0+(y.^3);
du_dy = 0+3.*x.*(y.^2);
d2u_dx2 = 0+0.*x;
d2u_dy2 = 0+6.*x.*y;
%}

%{
du_dt = 10+0.*x;
du_dx = 0+5.*x.^4;
du_dy = 0+5.*y.^4;
d2u_dx2 = 0+20.*x.^3;
d2u_dy2 = 0+20.*y.^3;
%}

%{
du_dt = 10+0.*x;
du_dx = 0+6.*x.^5;
du_dy = 0+6.*y.^5;
d2u_dx2 = 0+30.*x.^4;
d2u_dy2 = 0+30.*y.^4;
%}

lap = 0.*x;

%{
u_0 = 0.8; v_0 = 0.8;
alpha_x = 0.01; alpha_y = 0.01;
x_0 = 0.5; y_0 = 0.5;
%
a = 1/(1+4*time);

bb_1 = x-(u_0*time)-0.5;
bb_2 = (1+4*time)*alpha_x;
bb_3 = (bb_1).^2;
b = -bb_3/bb_2;

cc_1 = y-(v_0*time)-0.5;
cc_2 = (1+4*time)*alpha_y;
cc_3 = (cc_1).^2;
c = -cc_3/cc_2;
%
da_dt = -4/(1+4*time)^2;

trm_1 = 2*u_0*bb_1/bb_2;
trm_2 = 4*(bb_1).^2/(bb_2*(1+4*time));
db_dt = (trm_1+trm_2)*exp(b);

trm_3 = 2*v_0*cc_1/cc_2;
trm_4 = 4*(cc_1).^2/(cc_2*(1+4*time));
dc_dt = (trm_3+trm_4)*exp(c);
%
du_dt = (da_dt*exp(b)*exp(c))+(a*db_dt*c)+(a*b*dc_dt);
%
trm_1 = -2*bb_1/bb_2;
du_dx = a*exp(b)*exp(c)*trm_1;
%
trm_2 = -2*cc_1/cc_2;
du_dy = a*exp(b)*exp(c)*trm_2;
%
trm_ax = a*exp(c);
trm_bx_1 = -2/bb_2;
trm_bx_2 = 4*(bb_1.^2)/(alpha_x*bb_2);
trm_bx = (trm_bx_1+trm_bx_2)*exp(b);
%
d2u_dx2 = trm_ax*trm_bx;
%
trm_ay = a*exp(b);
trm_cy_1 = -2/cc_2;
trm_cy_2 = 4*(cc_1.^2)/(alpha_y*cc_2);
trm_cy = (trm_cy_1+trm_cy_2)*exp(c);
%
d2u_dy2 = trm_ay*trm_cy;

lap_1 = du_dt+(c_x.*du_dx)+(c_y.*du_dy)-nu.*(d2u_dx2+d2u_dy2); %the source term of the conv-diff eqn.
lap = 0.*lap_1;
%}

%disp('Hola')

%lap = 1;



