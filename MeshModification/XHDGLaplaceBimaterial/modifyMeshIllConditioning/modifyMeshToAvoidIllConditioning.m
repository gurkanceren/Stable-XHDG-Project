function [X,LS,movedNodes]=modifyMeshToAvoidIllConditioning(tol,Elements,X,T,LS,referenceElement)
% Element is assumed to be ill-conditioned if an area in the reference
% element is smaller than tol^2 (similar to distance of a node smaller than tol)

%Elements cut by the interface
Tcut = T(Elements.Int,:);
nOfElements = size(Tcut,1);

%Quadrature for standart triangle and quadrilateral
degree = referenceElement.degree;
[zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*degree);
wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
[zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*degree);

%nodal coordinates at reference element
nodesCoord = referenceElement.NodesCoord;
i=1;
movedNodes = []; aviso = 0;
%Loop to detect nodes close to the interface in ill-conditioned elements
while size(Tcut,1)>0
%for e=1:nOfElements
    %Te = Tcut(e,:);
    Te = Tcut(1,:);
    Xe = X(Te,:);
    LSe = LS(Te);
    
    [zgp,wgp,n1,n2,PtsInt] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
    A1 = sum(wgp(1:n1)); %area of region 1 (in reference element)
    A2 = sum(wgp(n1+1:end)); %area of region 2
    if (A1<tol^2 | A2<tol^2) %ill-conditioned element -> nodes moved to the interface
        disp('Ill-conditioned element detected:');
        Te
        LevelSet = LSe;
        movedNodesElem = [];
        nodesCloseInt = find(abs(LSe)<tol & LSe>0); %LS is assumed to be the signed distance
        for i=1:length(nodesCloseInt)
            inode = nodesCloseInt(i);
            Xi = nodesCoord(inode,:); %local coordinates
            if inode == 3 %third vertex (for triangles) has not proper definition of derivatives
                Xi = Xi+[1,-1]*tol*1.e-3;
            end
            Xp = normalProjectionTo0LS(Xi,Xe,LSe,referenceElement,tol); %projection with local coordinates
            if ~isempty(Xp)
                N = computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,Xp);
                N = N(:,:,1)';
                X(Te(inode),:) = N*Xe;
                LS(Te(inode)) = 0;
                movedNodesElem = [movedNodesElem,Te(inode)];
            end
            if isempty(movedNodesElem)
               aviso = 1; 
               disp('WARNING: an element is still ill-conditioned');
            else
                disp(sprintf('Fixed :-)\n\n'));
            end
            movedNodes = union(movedNodes,movedNodesElem);
        end
        if size(Tcut,1)>1 
            Tcut = Tcut(2:end,:);
            ElementsCut=SetElements(Tcut,LS,[1,0],referenceElement);
            Tcut = Tcut(ElementsCut.Int,:);
        else
            Tcut = [];
        end
    else
        if size(Tcut,1)>1, Tcut = Tcut(2:end,:); else, Tcut = []; end
    end
end

if aviso == 1
    disp('------------------------------------------------------------------------------')
    disp(sprintf('      [ WARNING: some ill-conditioned elements have not been fixed! ]')); 
    disp('------------------------------------------------------------------------------')
end

%__________________________________________
% Projection onto the 0-level set
function Xip = normalProjectionTo0LS(XiPoint,Xe,LSe,referenceElement,tol)
%Computes an approximation of the the normal projection of Xpoint to the
%0-level set
%The point is assumed to be close to the 0-level set
%The function returns [] if the proyection goes out from the element

N = computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,XiPoint);
Nxi = N(:,:,2)'; Neta = N(:,:,3)'; N = N(:,:,1)'; 

%Jacobian of the isoparametric transformation at this point (for 
%approximation of change of distance
J = [Nxi*Xe(:,1) Nxi*Xe(:,2);Neta*Xe(:,1) Neta*Xe(:,2)];

%value of level-set function and gradient at point Xpoint
LSval = N*LSe;
distance = LSval*sqrt(det(J)); %approximation of the distance to the 0-level set in local coordinates
LSgrad = [Nxi*LSe,Neta*LSe];
n = -sign(LSval)*LSgrad/norm(LSgrad);
P2 = XiPoint + 2*distance*n; %approximation to the symmetric point on the other side of the interface

N = computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,P2);
N = N(:,:,1)'; 
LSval2 = N*LSe;

Xip = FindPointLS(0,LSe,XiPoint,LSval,P2,LSval2,referenceElement);
%If the point is far from being inside the triangle the routine returns []
tolt = 0.1;
if (Xip(1) + Xip(2)>tolt || Xip(1)<-1-tolt || Xip(1)>1+tolt || Xip(2)<-1-tolt || Xip(2)>1+tolt)
    Xip = [];
end




