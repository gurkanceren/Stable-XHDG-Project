function L2Norm = computeL2Normpress(referenceElement,X,T,u,u0,time,Re,varargin)

% The function computeL2Norm allows to compute the L2 Norm using an
% analytical function handle. This function has to depend on x-y
% coordinates.
%
% Input:
%  referenceElement: information of the reference element
%  X: nodal coordinates
%  T: connectivity matrix
%  u: FEM solution (nodal)
%  u0: analytical function handle. This function has to have the x-y
%      coordinates matrix [xi yi] as a first input argument.
%  varargin: ordered list of input arguments needed to execute u0 but x-y
%            coordinates.
% Output:
%  L2Norm: computed L2 norm = sqrt(integral[(u-u0)^2]) over the domain defined by
%          T,X

%Number of elements and number of mesh nodes
nOfElements = size(T,1);


%Loop in 2D elements
L2Norm = 0;
for iElem = 1:nOfElements
    
    L2Norm = L2Norm + ElementalL2Normpress(X,referenceElement,u,u0,time,Re,varargin{:});
end
L2Norm = sqrt(L2Norm);


%_______________________________________________________________________
function elemL2Norm = ElementalL2Normpress(Xe,referenceElement,ue,u0,time,Re,varargin)

%Information of the reference element
IPw = referenceElement.IPweights;
N = referenceElement.N;
Nxi= referenceElement.Nxi;
Neta= referenceElement.Neta;
nOfElementNodes= size(N,2);

%Number of Gauss points
ngauss = length(IPw);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

%Compute elemental L2 Norm
elemL2Norm = 0;
for g = 1:ngauss
    %Values at current integration point
     N_g = N(g,:);
     Nxi_g = Nxi(g,:);
     Neta_g = Neta(g,:);
     
    xy_g = N_g*Xe;  % gauss pt coord
    
    ue_g = N_g*ue (1:nOfElementNodes);
    u0_g = u0(xy_g,time,Re,varargin{:});
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    %Integration weight
    dvolu=IPw(g)*det(J);
    %dvolu=IPw(g)*detJ;
    %Contribution of the current integration point to the elemental L2 Norm
    elemL2Norm = elemL2Norm + (ue_g-u0_g)^2*dvolu;
end






