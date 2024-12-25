function [u_star,shapeFunctions] = hdg_postprocess(X,T,u,q,referenceElement_star,referenceElement,Elements,mu)

nOfElements = size(T,1);
nOfElementNodes = size(T,2);
coordRef_star = referenceElement_star.NodesCoord;
npoints = size(coordRef_star,1);

% Vandermonde matrix
coordRef = referenceElement.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = referenceElement.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');

% Compute shape functions at interpolation points
shapeFunctions = zeros(npoints,nOfNodes);
for ipoint = 1:npoints
    p = orthopoly2D(coordRef_star(ipoint,:),nDeg);
    shapeFunctions(ipoint,:) = (invV*p)';
end

% u star initialization
u_star = zeros(npoints*nOfElements,1);

% Loop in elements
for iElem = 1:nOfElements
    
    Te_lin = T(iElem,:);
    Xold = X(Te_lin,:);
    Xe=shapeFunctions*Xold;
    LSe_star = EvaluateLS(Xe,1);
    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    ind_star = (iElem-1)*npoints+1:iElem*npoints;
    ue = shapeFunctions*u(ind);
    qe(1:2:2*npoints) = shapeFunctions*q(2*ind-1);
    qe(2:2:2*npoints) = shapeFunctions*q(2*ind);
    [Ke,Beqe,int_ue_star, int_ue] = ElementalMatrices(Xe,referenceElement_star,ue,qe,Elements,iElem,LSe_star,referenceElement,mu);
    
    % Lagrange multipliers
    K = [Ke int_ue_star; int_ue_star' 0];
    f = [Beqe;int_ue];
    
    % elemental solution
    sol = K\f;
    
    % postprocessed solution
    u_star(ind_star) = sol(1:end-1);
end

%%
%% ELEMENTAL MATRIX

function [K,Bq,int_u_star,int_u] = ElementalMatrices(Xe,referenceElement_star,ue,qe,Elements,iElem,LSe_star,referenceElement,mu)


c=find(iElem==Elements.D1);
d=find(iElem==Elements.D2);


nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);
K = zeros(nOfElementNodes_star,nOfElementNodes_star);
Bq = zeros(nOfElementNodes_star,1);
int_u_star = zeros(nOfElementNodes_star,1);
int_u = 0;


if ~isempty(c)
    
    % Information of the reference element
    IPw = referenceElement_star.IPweights;
    N = referenceElement_star.N;
    Nxi = referenceElement_star.Nxi;
    Neta = referenceElement_star.Neta;
    ngauss=size(IPw,1);
    weight=IPw;
    
    
elseif ~isempty(d)
    
    % Information of the reference element
    
    ngauss=0;
    IPw =0;
    N = 0;
    Nxi = 0;
    Neta = 0;
    NN_xy = 0 ;
    
else isempty(c) && isempty(d);
    
    
    p_star = referenceElement_star.degree;
    %Quadrature for standart triangle and quadrilateral
    [zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*p_star);
    wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
    [zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*p_star);
    
    
    [zgp,wgp,n1,n2,PtsInt,FaceInfo] = ModifyQuadrature(LSe_star,referenceElement_star,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement_star.degree,referenceElement_star.NodesCoord,zgp);
    N = shapeFunctions(:,:,1)';  Nxi = shapeFunctions(:,:,2)'; Neta= shapeFunctions(:,:,3)';
    weight=wgp; ngauss=n1;   %SFM
    
end


% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% VOLUME COMPUTATIONS: LOOP IN GAUSS POINTS
for g = 1:ngauss
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    
    %Integration weight
    dvolu=weight(g)*det(J);
    
    %x and y derivatives
    invJ = inv(J);
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    %     NN_xy(1:2:end) = Nx_g;
    %     NN_xy(2:2:end) = Ny_g;
    
    % u and q at gauss points
    u_g = N_g*ue;
    qx_g = N_g*qe(1:2:end)';
    qy_g = N_g*qe(2:2:end)';
    
    %Contribution of the current integration point to the elemental matrix
    Bq = Bq - (Nx_g'*qx_g + Ny_g'*qy_g)*dvolu/mu;
    K = K + (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    int_u_star = int_u_star + N_g'*dvolu;
    int_u = int_u + u_g*dvolu;
end


