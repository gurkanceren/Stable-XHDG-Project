function P0 = FindPointLS(linear,LSe,P1,d1,P2,d2,referenceElement)
%
% P0 = FindPointLS(linear,LSe,P1,d1,P2,d2,referenceElement)
% The function computes a point P0 on the interface and on the line defined by points P1 and P2 
% 
% Input:
%    linear: set this value to 1 to use a linear interpolation (usually not needed)
%    LSe: value of the level set function at the element's nodes
%    P1, P2: points defining the line where the interface point is searched. 
%    d1, d2: value of the level set function on P1 and P2
%    referenceElement


tol = 1e-20; 

nDeg = referenceElement.degree; 
Xe_ref = referenceElement.NodesCoord; 

v = P2 - P1;

if isempty(d1)
    d1 = FuncEvalLS(0,P1,v,LSe,nDeg,Xe_ref);
end
if isempty(d2)
    d2 = FuncEvalLS(0,P2,v,LSe,nDeg,Xe_ref);
end

if abs(d1) < tol 
    P0 = P1;
elseif abs(d2) < tol 
    P0 = P2;
elseif linear == 1  
    if d1*d2 > 0
        error('These two points are on the same subdomain')
    end
    
    lambda = d1/(d1-d2);
    P0 = P1 + lambda*v;
else
    if d1*d2 > 0
        lambda0 = 0; 
	else
        lambda0 = d1/(d1-d2);
    end
    
    options.TolX = tol;
    lambda = fzero(@(x)FuncEvalLS(x,P1,v,LSe,nDeg,Xe_ref),lambda0,options);
    P0 = P1 + lambda*v;    
    
    if (min(min(P0)) < -1) || (P0(2)+P0(1) > 0)
        lambda0 = 0.5; 
        lambda = fzero(@(x)FuncEvalLS(x,P1,v,LSe,nDeg,Xe_ref),lambda0,options);
        P0 = P1 + lambda*v;
    end
    
end




%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------


function resLS = FuncEvalLS(lambda,P,v,LSe,nDeg,Xe_ref)
%
% resLS is the value of the level-set function at the point P+lambda*v,
% in an element with nodal values of the level-set LSe.

N = computeShapeFunctionsAtPoints(nDeg,Xe_ref,P+lambda*v);
N = N(:,:,1)'; 
resLS = N*LSe;
