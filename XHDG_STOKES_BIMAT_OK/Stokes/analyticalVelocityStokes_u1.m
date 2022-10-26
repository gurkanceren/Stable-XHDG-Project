function u1=analyticalVelocityStokes_u1(X)

u=analyticalVelocityStokes(X);
n=size(u,2);
u1=u(:,1);