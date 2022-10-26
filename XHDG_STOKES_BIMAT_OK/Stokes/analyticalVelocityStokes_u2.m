function u2=analyticalVelocityStokes_u2(X)

u=analyticalVelocityStokes(X);
n=size(u,2);
u2=u(:,2);