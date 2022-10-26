function IntPt = intersectionPoint1D(LSe,referenceElement,tol)
% WARNING!: this only works for faces cut only at 1 point

%Checking if the element is cut at one point
if(LSe(1)*LSe(end)>0)
    error('This element is not cut at a unique point')
end

nDeg = referenceElement.degree;
Xi_ref = referenceElement.NodesCoord1d;

options.TolX = tol;
IntPt = fzero(@(x)FuncEvalLS(x,LSe,nDeg,Xi_ref),0,options);



%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

function resLS = FuncEvalLS(xi,LSe,nDeg,Xi_ref)
%
% resLS is the value of the level-set function at the point P+lambda*v,
% in an element with nodal values of the level-set LSe.

N = computeShapeFunctionsAtPoints(nDeg,Xi_ref,xi);
N = N(:,:,1)';
resLS = N*LSe;
