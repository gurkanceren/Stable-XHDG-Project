function [zgp,wgp,vect_mu,Pts] = ModifyQuadratureStdTri(Xe,LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua)
% 
% [zgp,wgp,vect_mu,Pts] = ModifyQuadratureStdTri(Xe,LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua)
% This function computes a quadrature rule for a triangle cut by the
% interface in an "standard way", that is, the interface cuts two edges,
% exactly one time each. 
% Appart from the integration points and weights, it returns a vector
% vect_mu indicating the subdomain to which each point velong and a vector
% Pts with p+1 points describing the interface. 
% Input: 
%    Xe: nodal coordinates
%    LSe: value of the level set function at the nodes
%    referenceElement
%    zgp_tri,wgp_tri,zgp_qua,wgp_qua: quadrature points and weights in a
%    reference triangle and quadrilateral.

tol = 1e-10;

nDeg = referenceElement.degree; 
Xe_ref = referenceElement.NodesCoord; 
AreaRefTri = sum(wgp_tri); 

mu1 =  1; 
mu2 = -1; 

LSe_Vertices = LSe(1:3);
if LSe_Vertices >= 0
    zgp = zgp_tri; wgp = wgp_tri; vect_mu = mu1*ones(length(wgp),1); Pts = []; 
elseif LSe_Vertices <= 0
    zgp = zgp_tri; wgp = wgp_tri; vect_mu = mu2*ones(length(wgp),1);  Pts = []; 
else
    Vert0 = find(abs(LSe_Vertices)<tol);
    if length(Vert0) == 1
        Pts = FindPointsInterface(LSe, referenceElement);
        node0 = Vert0; 
        node1 = node0+1; if(node1) == 4, node1 = 1; end
        node2 = node0-1; if(node2) == 0, node2 = 3; end
        if LSe(node1,:) > 0
            mu_1 = mu1; mu_2 = mu2; 
        else
            mu_1 = mu2; mu_2 = mu1; 
        end
        [zgp1,wgp1,vect_mu1] = ModQuadTri(Pts,Xe_ref(node1,:),zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu_1,referenceElement);
        [zgp2,wgp2,vect_mu2] = ModQuadTri(Pts,Xe_ref(node2,:),zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu_2,referenceElement);
        zgp = [zgp1; zgp2]; 
        wgp = [wgp1, wgp2]; 
        vect_mu = [vect_mu1; vect_mu2]; 
    else
        aux = find(LSe_Vertices > 0);
        n1 = length(aux); 
        if n1 == 1
            mu_1 = mu1; mu_2 = mu2; 
            node0 = find(LSe_Vertices > 0);
            node1 = node0+1; if node1 == 4, node1 = 1; end
            node2 = node0-1; if node2 == 0, node2 = 3; end
        else
            mu_1 = mu2; mu_2 = mu1; 
            node0 = find(LSe_Vertices <= 0);
            node1 = node0+1; if node1 == 4, node1 = 1; end
            node2 = node0-1; if node2 == 0, node2 = 3; end
        end
        
        Pts = FindPointsInterface(LSe, referenceElement);
        
        [zgp1,wgp1,vect_mu1] = ModQuadTri(Pts,Xe_ref(node0,:),zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu_1,referenceElement);
        Pts_aux = [Xe_ref(node2,:); Xe_ref(node1,:)]; 
        [zgp2,wgp2,vect_mu2] = ModQuadQua(Pts,Pts_aux,zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu_2,referenceElement);

        zgp = [zgp1; zgp2]; 
        wgp = [wgp1, wgp2]; 
        vect_mu = [vect_mu1; vect_mu2];         
    end
end


if ~isempty(zgp)
    N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,Pts);
    N = N(:,:,1)'; 
    Pts = N*Xe;
    Pts = reshape(Pts', 1,2*size(Pts,1)); 

    N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,zgp);
    N = N(:,:,1)'; 
    zgp = N*Xe;

    T = Xe(1:3,:); 
    v1 = [T(1,:)-T(3,:),0];
    v2 = [T(2,:)-T(3,:),0];
    AreaTri = norm(cross(v1,v2),2)/2;
    wgp = wgp*AreaTri/AreaRefTri;
end
