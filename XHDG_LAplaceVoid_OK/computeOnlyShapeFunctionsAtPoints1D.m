%____________________________________________________________________
function shapeFunctions = computeOnlyShapeFunctionsAtPoints1D(z,nDeg,coord)
%number of nodes/polynomials
nOfNodes = nDeg+1;
if nOfNodes~=size(coord,1)
    error('Error locatePointscurvedSegmentElement')
end

nOfPoints = length(z);

%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

shapeFunctions = zeros(nOfPoints,nOfNodes);

%Integration over [-1,1]
for i = 1:nOfPoints
    x = z(i);
    p = orthopoly1D(x,nDeg);
    N = U\(L\(P*p));
    shapeFunctions(i,:) = N';
end

