function [zgp,wgp,vect_mu] = ModQuadQua(Pts_int,Pts_aux,zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu,referenceElement)
% 
% [zgp,wgp,vect_mu] = ModQuadQua(Pts_int,Pts_aux,p,zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu,referenceElement)
% 
% This function computes a quadrature rule for a quadrilateral cell with a curved edge. 
% Input: 
%    Pts_int: describe the curved edge
%    Pts_aux: are the other two vertices in the quadrilateral
%    zgp_tri,wgp_tri,zgp_qua,wgp_qua: quadrature rules in a reference element
%    mu: material properties
%    referenceElement


tol = norm(Pts_int(1,:) - Pts_int(end,:))*1e-20;

p = referenceElement.degree; 
Xe_ref_1D = referenceElement.NodesCoord1d; 

if norm(Pts_int(1,:)-Pts_aux(1,:))<tol ||  norm(Pts_int(end,:)-Pts_aux(1,:))<tol
    [zgp,wgp,vect_mu] = ModQuadTri(Pts_int,Pts_aux(:,2),zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu,referenceElement); 
elseif norm(Pts_int(1,:)-Pts_aux(2,:))<tol ||  norm(Pts_int(end,:)-Pts_aux(2,:))<tol
    [zgp,wgp,vect_mu] = ModQuadTri(Pts_int,Pts_aux(:,1),zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu,referenceElement); 
else
    [Pts_int,Pts_aux] = CheckOrientation(Pts_int,Pts_aux);

    Pts_aux = SetPtsOnLine(Pts_aux([2,1],:),p,Xe_ref_1D);
    
    [aux, ind] = sort(zgp_qua(:,1)); 
    zgp_qua = zgp_qua(ind,:); 
    wgp_qua = wgp_qua(ind); 
    r = zgp_qua(:,1); s = zgp_qua(:,2); 
    
    NT = computeShapeFunctionsAtPoints(p,Xe_ref_1D,r); 
    N   = NT(:,:,1)'; 
    Nxi = NT(:,:,2)'; 
    
    x = 0.5*(  (1+s).*(N*Pts_aux(:,1)) + (1-s).*(N*Pts_int(:,1))  ); 
    y = 0.5*(  (1+s).*(N*Pts_aux(:,2)) + (1-s).*(N*Pts_int(:,2))  ); 
    zgp = [x,y]; 

    ngaus = length(wgp_qua);
    wgp = zeros(1,ngaus); vect_mu = mu*ones(ngaus,1);

    dx_r = 0.5*(  (1+s).*(Nxi*Pts_aux(:,1)) + (1-s).*(Nxi*Pts_int(:,1))  ); 
    dx_s = 0.5*(  (N*Pts_aux(:,1)) - (N*Pts_int(:,1))  ); 
    dy_r = 0.5*(  (1+s).*(Nxi*Pts_aux(:,2)) + (1-s).*(Nxi*Pts_int(:,2))  ); 
    dy_s = 0.5*(  (N*Pts_aux(:,2)) - (N*Pts_int(:,2))  );     
    for i = 1:ngaus
        J = [dx_r(i)   dy_r(i);  dx_s(i)   dy_s(i)];
        wgp(i) = wgp_qua(i)*det(J); 
    end
end



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



function  [Pts_int,Pts_aux] = CheckOrientation(Pts_int,Pts_aux)

P1 = Pts_int(1,:);
P2 = Pts_int(end,:);
P3 = Pts_aux(1,:);
P4 = Pts_aux(2,:); 

v1 = P2-P1; v1 = v1/norm(v1);
v2 = P3-P1; v2 = v2/norm(v2);
v3 = P4-P1; v3 = v3/norm(v3);

v_1 = cross([v1,0], [v2,0]);
v_2 = cross([v2,0], [v3,0]);
if v_1(3) < 0
    Pts_int = [flipud(Pts_int(:,1)), flipud(Pts_int(:,2))];
    if v_2(3) < 0
        Pts_aux = [flipud(Pts_aux(:,1)), flipud(Pts_aux(:,2))];
    end
end




function Pts = SetPtsOnLine(Pts,p,lambda)

if size(Pts,1) == 2 && p > 1
    P1 = Pts(1,:);
    P2 = Pts(2,:);
    v = (P2-P1); 
    
    lambda = (lambda - lambda(1)) / (lambda(end) - lambda(1)); 
    Pts = [P1(1)+lambda*v(1), P1(2)+lambda*v(2)];
end