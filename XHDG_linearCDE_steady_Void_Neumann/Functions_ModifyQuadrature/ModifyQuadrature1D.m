function [zgp,wgp,n1,n2,IntPt] = ModifyQuadrature1D(LSe,referenceElement)
% [zgp,wgp,n1,n2,PtsInt_e] = ModifyQuadrature1D(LSe,referenceElement)
% Input: LS values at the element, reference element, standard 1D
% quadrature
% Output: quadrature for the cut element, number of integration points for
% each domain, and Intersection point


IntPt = intersectionPoint1D(LSe,referenceElement,1.e-16);

w = referenceElement.IPweights1d;
z = referenceElement.IPcoordinates1d;

vertex1 = referenceElement.NodesCoord1d(1);
vertex2 = referenceElement.NodesCoord1d(end);

[z1,w1]=segmentQuadrature1D(vertex1,IntPt,z,w);
[z2,w2]=segmentQuadrature1D(IntPt,vertex2,z,w);

if LSe(1)>0 %1st segment is domain 1
    n1 = length(w1); n2 = length(w2);
    zgp = [z1;z2]; wgp = [w1;w2];
else %1st segment is domain 2
    n2 = length(w1); n1 = length(w2);
    zgp = [z2;z1]; wgp = [w2;w1];
end


%______________________________________________________
% Given de points a and b at the reference element, returns a quadrature in
% the segment [a,b]
function [zs,ws]=segmentQuadrature1D(a,b,z,w)

l = (b-a)/2;
ws = w*l;
zs = (a+b)/2 + z*l;
