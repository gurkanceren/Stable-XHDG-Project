function gi = bcNeumannPcorr(Nx,Ny,x,y)

%{
    Re = 10;
    a = Re*Re/4;
    b = 4*pi*pi;
    lmbda = (Re/2)-sqrt(a+b);

    du_dx = -2*lmbda*exp(2*lmbda*x).*cos(pi*(4*y-1));
    du_dy = 4*pi*exp(2*lmbda*x).*sin(pi*(4*y-1));
    dv_dx = (lmbda*lmbda/pi)*exp(2*lmbda*x).*sin(pi*(4*y-1));
    dv_dy = 2*lmbda*exp(2*lmbda*x).*cos(pi*(4*y-1));
    press = (-1/2)*exp(4*lmbda*x);
%}

%
    du_dx = 0.*x;
    du_dy = 0.*y;
    dv_dx = 0.*x;
    dv_dy = 0.*y;
    press = 0.*(x+y);
%
    gx = (-mu*du_dx.*Nx)-(mu*dv_dx.*Ny)+(press.*Nx);
    gy = (-mu*du_dy.*Nx)-(mu*dv_dy.*Ny)+(press.*Ny);

    gi = [gx , gy]; 

