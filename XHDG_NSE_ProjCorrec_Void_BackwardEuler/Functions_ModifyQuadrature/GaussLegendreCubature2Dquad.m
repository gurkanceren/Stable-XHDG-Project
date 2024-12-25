function [xi,w_g] = GaussLegendreCubature2Dquad(order)
%
% [xi,w_g] = GaussLegendreCubature2Dquad(order)
%
% Input:    order:        order for the cubature rule 
% Output:   xi, w_g:      integration points and weights in teh reference element



n = ceil((order+1)/2); 
ngaus = n+1;
[pg_1D, wg_1D] = gaussLegendre(ngaus,-1,1);

xx = pg_1D';
yy = pg_1D';
[xx,yy] = meshgrid(xx,yy);
xx = reshape(xx,ngaus^2,1);
yy = reshape(yy,ngaus^2,1);

xi = [xx,yy];
w_g = wg_1D'*wg_1D;
w_g = reshape(w_g,1,ngaus^2);

    