function u1=analyticalVelocityStokes_u1(X,time,Re)

u=analyticalVelocityStokes(X,time,Re);
n=size(u,1);
u1=u(1:n/2);