function [zgp,wgp,vect_mu] = ModQuadTri(Pts_int,P0,zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu,referenceElement)
% 
% [zgp,wgp,vect_mu] = ModQuadTri(Pts_int,P0,p,zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu,referenceElement)
% 
% This function computes a quadrature rule for a triangular cell with a curved edge. 
% Input: 
%    Pts_int: describe the curved edge
%    P0: is the vertex to complete the triangle
%    zgp_tri,wgp_tri,zgp_qua,wgp_qua: quadrature rules in a reference element
%    mu: material properties
%    referenceElement


tol = norm(Pts_int(1,:) - Pts_int(end,:))*1e-10;

p = referenceElement.degree; 
Xe_ref_1D = referenceElement.NodesCoord1d; 

if size(Pts_int(1,:),1)==size(P0,1)
    P0=P0;
else
    P0=P0';
end

if norm(Pts_int(1,:)-P0) < tol || norm(Pts_int(end,:)-P0) < tol
    zgp = []; wgp = []; vect_mu = [];
else
    [Pts_int,P0] = CheckOrientation(Pts_int,P0);

    Tri_aux = [Pts_int(1,:); Pts_int(end,:); P0];
    v1 = [Tri_aux(1,:)-Tri_aux(3,:),0];
    v2 = [Tri_aux(2,:)-Tri_aux(3,:),0];
    Area = norm(cross(v1,v2),2)/2;
    
    if abs(Area) < 1e-20
        zgp = []; wgp = []; vect_mu = [];
    else
        Pts_int = SetPtsOnLine(Pts_int,p,Xe_ref_1D);

        r = zgp_qua(:,1); s = zgp_qua(:,2);
        NT = computeShapeFunctionsAtPoints(p,Xe_ref_1D,r); 
        N   = NT(:,:,1)'; 
        Nxi = NT(:,:,2)'; 

        x = 0.5*((1+s).*P0(1) + (1-s).*(N*Pts_int(:,1)));
        y = 0.5*((1+s).*P0(2) + (1-s).*(N*Pts_int(:,2)));
        zgp = [x,y];

        ngaus = length(wgp_qua);
        wgp = zeros(1,ngaus); vect_mu = mu*ones(ngaus,1);

        dx_r = 0.5*((1-s).*(Nxi*Pts_int(:,1)));
        dx_s = 0.5*(P0(1) - (N*Pts_int(:,1)));
        dy_r = 0.5*((1-s).*(Nxi*Pts_int(:,2)));
        dy_s = 0.5*(P0(2) - (N*Pts_int(:,2)));
        for i = 1:ngaus
            J = [dx_r(i)   dy_r(i);  dx_s(i)   dy_s(i)];
            wgp(i) = wgp_qua(i)*det(J);
        end
    end
end



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------




function [Pts_int,P0] = CheckOrientation(Pts_int,P0)

P1 = Pts_int(1,:);
P2 = Pts_int(end,:);
P3 = P0;

v1 = P2-P1; v1 = v1/norm(v1);
v2 = P3-P1; v2 = v2/norm(v2);

v = cross([v1,0], [v2,0]);
if v(3) < 0
    Pts_int = [flipud(Pts_int(:,1)), flipud(Pts_int(:,2))];
end




function Pts = SetPtsOnLine(Pts,p,lambda)

if size(Pts,1) == 2 && p > 1
    P1 = Pts(1,:);
    P2 = Pts(2,:);
    v = (P2-P1); 
    
    lambda = (lambda - lambda(1)) / (lambda(end) - lambda(1)); 
    Pts = [P1(1)+lambda*v(1), P1(2)+lambda*v(2)];
end
