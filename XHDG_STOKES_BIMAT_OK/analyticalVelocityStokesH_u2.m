function u2=analyticalVelocityStokesH_u2(X)

u=analyticalVelocityStokesH(X);
n=size(u,2);
u2=u((2*n/4)+1:3*n/4);