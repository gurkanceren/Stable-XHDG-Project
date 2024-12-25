function L2Norm = computeL2Norm_cut(LSe,referenceElement,X,T,u,u0,time,Re,varargin)

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

%Shape Functions gauss points and weights for cut elements
p = referenceElement.degree;
%Quadrature for standart triangle and quadrilateral
[zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*p);
wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
[zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*p);


[zgp,wgp,n1,n2,PtsInt,FaceInfo] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zgp);
N = shapeFunctions(:,:,1)';  Nxi = shapeFunctions(:,:,2)'; Neta= shapeFunctions(:,:,3)';
weight=wgp; ngauss=n1;


%Loop in 2D elements
L2Norm = 0;
for iElem = 1:nOfElements
    L2Norm = L2Norm + ElementalL2Norm(X,u,u0,N,Nxi,Neta,weight,ngauss,time,Re,varargin{:});
end
L2Norm = sqrt(L2Norm);


%_______________________________________________________________________
function elemL2Norm = ElementalL2Norm(Xe,u,u0,N,Nxi,Neta,weight,ngauss,time,Re,varargin)

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);


%Compute elemental L2 Norm
elemL2Norm = 0;
for g = 1:ngauss
    %Values at current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    xy_g = N_g*Xe;   
    ue_g = N_g*u;
    u0_g = u0(xy_g,time,Re,varargin{:});
   
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe     Neta_g*ye];
    %Integration weight
    dvolu=weight(g)*det(J);
    %Contribution of the current integration point to the elemental L2 Norm
    elemL2Norm = elemL2Norm + (ue_g-u0_g)^2*dvolu;
    
end





