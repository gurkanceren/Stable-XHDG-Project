function [zgp,wgp,n1,n2,PtsInt_e,CutFaces_e] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua)
%
% [zgp,wgp,n1,n2,PtsInt_e] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua)
% Input: 
%     LSe: value of the level set function at the element's nodes
%     referenceElement: properties of the element used
%     zgp_tri, wgp_tri: quadrature in the reference triangle
%     zgp_qua, wgp_qua: quadrature in the reference quadrilateral
% Output: 
%     zgp,wgp: numerical quadrature in the element cut by the interface
%     n1,n2: number of integration points in each subdomain. 
%            Points are reordered so that zgp(1:n1,:) are the integration
%            points in subdomain 1 while zgp(n1+(1:n2),:) are the ones in
%            subdomain 2. 
%     PtsInt_e: list of points defining the interface in the reference element
%               It is a matrix with 2(p+1) columns and k rows. 
%               Each row is a set of p+1 points defining a p-th degree
%               polynimial that approximates a piece of interface: 
%               [x1,y1, x2,y2, ... x_{p+1},y_{p+1}]



nen = length(LSe);
n = 2*nen;

% If the element is cut in an standard way (the interface cuts two edges
% exactly once) we computed the modified quadrature. 
% If not, the element is divided.
if CheckStdElement(LSe,referenceElement)
    [zgp,wgp,vect_mu,PtsInt_e,CutFaces_e] = SecondVersion_of_ModifyQuadratureStdTri(referenceElement.NodesCoord,LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
else
    warning('Element cut faces are not correctly defined')
    CutFaces_e = [];
    [Xe1,LSe1,Xe2,LSe2] = SplitTri2Tri(referenceElement.NodesCoord,LSe,referenceElement);
    List_NoStdTri = [reshape(Xe1',1,n), LSe1'; reshape(Xe2',1,n), LSe2'];
    
    zgp = []; wgp = []; vect_mu = []; PtsInt_e = [];
    while ~isempty(List_NoStdTri)
        Xe_aux = reshape(List_NoStdTri(1,1:n),2,nen)';
        LSe_aux = List_NoStdTri(1,n+1:end)';
        if CheckStdElement(LSe_aux,referenceElement)
            [zgp1,wgp1,vect_mu1,PtsInt_e1,CutFaces_e1] = SecondVersion_of_ModifyQuadratureStdTri(Xe_aux,LSe_aux,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
            zgp = [zgp; zgp1];
            wgp = [wgp, wgp1];
            vect_mu = [vect_mu; vect_mu1];
            PtsInt_e = [PtsInt_e; PtsInt_e1];
            CutFaces_e = [CutFaces_e; CutFaces_e1]; 
        else
            [Xe1,LSe1,Xe2,LSe2] = SplitTri2Tri(Xe_aux,LSe,referenceElement);
            List_NoStdTri = [List_NoStdTri; ...
                reshape(Xe1',1,n), LSe1'; reshape(Xe2',1,n), LSe2'];
        end
        List_NoStdTri(1,:) = [];
    end
end

ind1 = find(vect_mu ==  1); n1 = length(ind1);
ind2 = find(vect_mu == -1); n2 = length(ind2);
zgp = zgp([ind1;ind2],:);
wgp = wgp([ind1;ind2]);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


function [Xe1,LSe1,Xe2,LSe2] = SplitTri2Tri(Xe,LSe,referenceElement)
%
% [Xe1,LSe1,Xe2,LSe2] = SplitTri2Tri(Xe,LSe,referenceElement)


P1 = Xe(1,:); P2 = Xe(2,:); P3 = Xe(3,:);

Vect_CutPtsEdge  = CountCutPoints(LSe,referenceElement);
nCutPoints = sum(Vect_CutPtsEdge);
nCutEdges = length(find(Vect_CutPtsEdge ~= 0));
if nCutPoints == 2 && nCutEdges == 1
    ind = find(Vect_CutPtsEdge ~= 0);
else
    h = [norm(P2-P1), norm(P3-P2), norm(P3-P1)];
    [aux,ind] = max(h);
end

if ind == 1
    PM = (P1+P2)/2;
    P1_1 = P2; P2_1 = P3;
    P1_2 = P3; P2_2 = P1;
elseif ind == 2
    PM = (P2+P3)/2;
    P1_1 = P3; P2_1 = P1;
    P1_2 = P1; P2_2 = P2;
else
    PM = (P3+P1)/2;
    P1_1 = P1; P2_1 = P2;
    P1_2 = P2; P2_2 = P3;
end

nDeg = referenceElement.degree;
Xe_ref = referenceElement.NodesCoord;
Xe_ref_lin = referenceElement.NodesCoord(1:3,:);
N = computeShapeFunctionsAtPoints(1,Xe_ref_lin,Xe_ref);
N = N(:,:,1)';
Xe1 = N*[PM; P1_1; P2_1];
Xe2 = N*[PM; P1_2; P2_2];

N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,Xe1); N = N(:,:,1)'; LSe1 = N*LSe;
N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,Xe2); N = N(:,:,1)'; LSe2 = N*LSe;

