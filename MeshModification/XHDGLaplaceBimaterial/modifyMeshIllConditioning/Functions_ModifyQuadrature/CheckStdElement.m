function res = CheckStdElement(LSe,referenceElement)
% 
% res = CheckStdElement(LSe,referenceElement)
% The function only works for triangular elements. 
% Input: 
%     LSe: nodal values of the level set function
%     referenceElement
% Output: 
%     res = 1 if the element is standard (that is, the interface cuts 2
%     edges, exactly once each)
%     res = 0 otherwise

tol = 1e-20; 

Vect_CutPtsEdge = CountCutPoints(LSe,referenceElement);
nCutPoints = sum(Vect_CutPtsEdge);
nCutEdges = length(find(Vect_CutPtsEdge ~= 0));

nen = length(LSe);
nDeg = referenceElement.degree; 
Xe_ref = referenceElement.NodesCoord; 


npt = 5*nDeg+1;
x = linspace(-1,1,npt);
n = npt*(npt+1)/2;
pts = zeros(n,2);
ini = 0;
for i = 1:npt
    ind = ini + (1:npt-i+1);
    pts(ind,1) = x(1:npt-i+1)';
    pts(ind,2) = x(i);
    ini = i*(npt - 0.5*(i-1));
end
N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,pts);
N = N(:,:,1)'; 
LS_aux = N*LSe;
aux = find(abs(LS_aux)<tol);
LS_aux(aux) = zeros(size(aux));

if all(LS_aux >= 0) || all(LS_aux <= 0) || (nCutEdges == 2 && nCutPoints == 2)
    res = 1;
else
    res = 0;
end



