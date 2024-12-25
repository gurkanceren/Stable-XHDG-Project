function u1=analyticalVelocityStokes_u1(X,time)

u=analyticalVelocityStokes(X,time);
n=size(u,1);
u1=u(1:n/2);