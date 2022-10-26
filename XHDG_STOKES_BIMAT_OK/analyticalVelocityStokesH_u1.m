function u1=analyticalVelocityStokesH_u1(X)

u=analyticalVelocityStokesH(X);
n=size(u,2);
u1=u(1:n/4);