function [zgp,wgp,vect_mu,Pts,CutFaces] = SecondVersion_of_ModifyQuadratureStdTri(Xe,LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua)
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

tol = 1e-20;

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
        %plot(Pts(:,1),Pts(:,2),'*r');
        Pts=locatePointscurvedSegmentElement(referenceElement.NodesCoord1d,Pts);
        %plot(Pts(:,1),Pts(:,2),'oc');
        
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
        figure(1); hold on;
        
        %plot(Pts(:,1),Pts(:,2),'*r');
        
        Pts=locatePointscurvedSegmentElement(referenceElement.NodesCoord1d,Pts);
        %plot(Pts(:,1),Pts(:,2),'oc');
        
        [zgp1,wgp1,vect_mu1] = ModQuadTri(Pts,Xe_ref(node0,:),zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu_1,referenceElement);
        Pts_aux = [Xe_ref(node2,:); Xe_ref(node1,:)]; 
        [zgp2,wgp2,vect_mu2] = ModQuadQua(Pts,Pts_aux,zgp_tri,wgp_tri,zgp_qua,wgp_qua,mu_2,referenceElement);

        zgp = [zgp1; zgp2]; 
        wgp = [wgp1, wgp2]; 
        vect_mu = [vect_mu1; vect_mu2];         
    end
end

CutFaces = [];
if ~isempty(zgp)
     if ~isempty(Pts)
    Pts = CorrectOrientation(Pts,LSe,referenceElement);    
    N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,Pts);
    N = N(:,:,1)'; 
    Pts_aux = N*Xe;
    Pts = reshape(Pts_aux', 1 ,(2*(nDeg+1))); 
    CutFaces = FindCutFaces(Pts_aux, referenceElement); 
    end
    
         
    
    N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,zgp);
    N = N(:,:,1)'; 
    zgp = N*Xe;

    T = Xe(1:3,:); 
    v1 = [T(1,:)-T(3,:),0];
    v2 = [T(2,:)-T(3,:),0];
    AreaTri = norm(cross(v1,v2),2)/2;
    wgp = wgp*AreaTri/AreaRefTri;
end



function CutFaces = FindCutFaces(Pts, referenceElement)

tol = 1e-10; 

Xe_ref = referenceElement.NodesCoord; 
faceNodes = referenceElement.faceNodes; 
n = size(faceNodes,2); 
PointsFace1 = Xe_ref(faceNodes(1,[1,n]),:); v1 = PointsFace1(2,:) - PointsFace1(1,:); 
PointsFace2 = Xe_ref(faceNodes(2,[1,n]),:); v2 = PointsFace2(2,:) - PointsFace2(1,:); 
PointsFace3 = Xe_ref(faceNodes(3,[1,n]),:); v3 = PointsFace3(2,:) - PointsFace3(1,:); 
Pref1 = Xe_ref(1,:); Pref2 = Xe_ref(2,:); Pref3 = Xe_ref(3,:); 

P1 = Pts(1,:); 
if norm(Pref1 - P1) < tol || norm(Pref2 - P1) < tol || norm(Pref3 - P1) < tol
    face1 = []; 
else
    M1 = [v1; P1 - PointsFace1(1,:)]; 
    if det(M1) < tol
        face1 = 1; 
    else
        M2 = [v2; P1 - PointsFace2(1,:)]; 
        if det(M2) < tol
            face1 = 2; 
        else
            M3 = [v3; P1 - PointsFace3(1,:)]; 
            if det(M3) < tol
                face1 = 3; 
            else
                face1 = 0; 
            end
        end
    end
end

P2 = Pts(end,:); 
if norm(Pref1 - P2) < tol || norm(Pref2 - P2) < tol || norm(Pref3 - P2) < tol
    face2 = []; 
else
    M1 = [v1; P2 - PointsFace1(1,:)]; 
    if det(M1) < tol
        face2 = 1; 
    else
        M2 = [v2; P2 - PointsFace2(1,:)]; 
        if det(M2) < tol
            face2 = 2; 
        else
            M3 = [v3; P2 - PointsFace3(1,:)]; 
            if det(M3) < tol
                face2 = 3; 
            else
                face2 = 0; 
            end
        end
    end
end
if isempty(face1)
        CutFaces = [face2, P2]; 
elseif isempty(face2)
    CutFaces = [face1, P1]; 
else
    CutFaces = [face1, P1; face2, P2]; 
end




function Pts = CorrectOrientation(Pts,LSe,referenceElement)
p = referenceElement.degree;
Xe_ref = referenceElement.NodesCoord;
n = floor((p+1)/2);
P1 = Pts(n,:);
P2 = Pts(n+1,:);
dx = P2-P1;
nt = [dx(2),-dx(1)];
P = P1 + nt/10;
N = computeShapeFunctionsAtPoints(p,Xe_ref,P);
N = N(:,:,1)';
d = N*LSe;
if d > 0
Pts = [flipud(Pts(:,1)), flipud(Pts(:,2))];
end

