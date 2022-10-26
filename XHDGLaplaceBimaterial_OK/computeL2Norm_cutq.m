function L2Norm = computeL2Norm_cutq(LSe,referenceElement,X,T,q,varargin)

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
Nold = shapeFunctions(:,:,1)';  Nxi_old = shapeFunctions(:,:,2)'; Neta_old= shapeFunctions(:,:,3)';
weight=wgp; ngauss=n1+n2;


%Enrichment of Shape Functions 
Nl=Nold(1:n1,:)*-1;
Nr=Nold(n1+1:n2+n1,:)*1;
NH=[Nl;Nr];
N=[Nold,NH];

Nlxi=Nxi_old(1:n1,:)*-1;
Nrxi=Nxi_old(n1+1:n2+n1,:)*1;
NHxi=[Nlxi;Nrxi];
Nxi=[Nxi_old,NHxi];

Nleta=Neta_old(1:n1,:)*-1;
Nreta=Neta_old(n1+1:n2+n1,:)*1;
NHeta=[Nleta;Nreta];
Neta=[Neta_old,NHeta];

%Loop in 2D elements
L2Norm = 0;
for iElem = 1:nOfElements
    L2Norm = L2Norm + ElementalL2Norm(X,q,N,Nxi,Neta,Nold,Nxi_old,Neta_old,weight,ngauss,varargin{:});
end
L2Norm = sqrt(L2Norm);


%_______________________________________________________________________
function elemL2Norm = ElementalL2Norm(Xe,qe,N,Nxi,Neta,Nold,Nxi_old,Neta_old,weight,ngauss,varargin)

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);


%Compute elemental L2 Norm
elemL2Norm = 0;
for g = 1:ngauss
    %Values at current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    N_g_old = Nold(g,:);
    Nxi_g_old = Nxi_old(g,:);
    Neta_g_old = Neta_old(g,:);
    
    xy_g = N_g_old*Xe;   
    qx_g =N_g*qe(1:2:end);
    %qy_g =N_g*qe(2:2:end);
    [q0_gx, q0_gy]=analiticalSolutionq(xy_g);
   
    %Jacobian
    J = [Nxi_g_old*xe	  Nxi_g_old*ye
        Neta_g_old*xe     Neta_g_old*ye];
    %Integration weight
    dvolu=weight(g)*det(J);
    %Contribution of the current integration point to the elemental L2 Norm
     elemL2Norm = elemL2Norm + (qx_g-q0_gx)^2*dvolu;
end






