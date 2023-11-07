function Xp = normalProjectionTo0LS(Xpoint,LSe,referenceElement,tol)
%Computes an approximation of the the normal projection of Xpoint to the
%0-level set
%The point is assumed to be close to the 0-level set, with a distance in
%the reference element smaller than tol

N = computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,Xpoint);
Nxi = N(:,:,2)'; Neta = N(:,:,3)'; N = N(:,:,1)'; 

%value of level-set function and gradient at point Xpoint
LSval = N*LSe;
LSgrad = [Nxi*LSe,Neta*LSe];
n = -sign(LSval)*LSgrad/norm(LSgrad);
P2 = Xpoint + 2*tol*n; %point symmetric to Xpoint on the other side of the interface

N = computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,P2);
N = N(:,:,1)'; 
LSval2 = N*LSe;

Xp = FindPointLS(0,LSe,Xpoint,LSval,P2,LSval2,referenceElement);





