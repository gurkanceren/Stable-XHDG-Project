function N=computeShapeFunctionsAtPoints(nDeg,coord,points,varargin)
%
% N=computeShapeFunctionsReferenceElement(nDeg,coord,points,varargin)
%
% Function to compute the shape functions (& derivatives) at Gauss points
%
% Input:
% nDeg:  degree of interpolation
% coord: nodal coordinates at the reference element
% points: points for evaluation of shape functions
% elementType (optional): 0 for quadrilateral, 1 for triangle. If it isn't
%                         given only triangle or 1D elements are considered
%
% Output:
% N: shape functions evaluated at the points (npoints x nnodes)
%

if ~isempty(varargin)
    elementType = varargin{:};
    if elementType == 0
        shapeFunctions=computeShapeFunctionsAtPointsQua(nDeg,coord,points);
        N = shapeFunctions(:,:,1)';
        return
    end
end

nsd = size(coord,2);
if nsd==1
    N=computeShapeFunctionsPoints1D(nDeg,coord,points);
elseif nsd==2
    N=computeShapeFunctionsPoints2D(nDeg,coord,points);
elseif nsd==3
    N=computeShapeFunctionsPoints3D(nDeg,coord,points);
else
    error('wrong nsd in computeShapeFunctionsPoints')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=computeShapeFunctionsPoints1D(nDeg,coord,points)

%number of nodes/polynomials
nOfNodes = nDeg+1;
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsPints1D')
end

z = points; nOfPoints = size(z,1);

%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

N = zeros(nOfPoints,nOfNodes);

%Integration over [-1,1]
for i = 1:nOfPoints
    x = z(i);
    p = orthopoly1D(x,nDeg);
    N(i,:)=U\(L\(P*p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=computeShapeFunctionsPoints2D(nDeg,coord,points)

%number of nodes/polynomials

nOfNodes = (nDeg+1)*(nDeg+2)/2;
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsPoints2D')
end

z = mapXiEtaZeta2rst(points);
nOfPoints = size(z,1);
%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

N = zeros(nOfPoints,nOfNodes);

for i = 1:nOfPoints
    x = z(i,:); %(r,s) coordinates
    p = orthopoly2D_rst(x,nDeg);
    N(i,:) = U\(L\(P*p));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=computeShapeFunctionsPoints3D(nDeg,coord,points)

%number of nodes/polynomials
nOfNodes = (nDeg+1)*(nDeg+2)*(nDeg+3)/6;
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsPoints2D')
end

z = mapXiEtaZeta2rst(points);
nOfPoints = size(z,1);
%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

N = zeros(nOfPoints,nOfNodes);

for i = 1:nOfPoints
    x = z(i,:); % (r,s,t) coordinates
    p = orthopoly3D_rst(x,nDeg);
    N(i,:) = U\(L\(P*p));
end

