function u2=analyticalVelocityStokes_u2_old(X,time)

u=analyticalVelocityStokes(X,time);
n=size(u,1);
u2=u(n/2+1:end);