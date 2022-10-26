
function [u_star,N] = hdg_postprocess4(X,T,u,q,referenceElement_star,referenceElement,Elements,mu1,mu2)

nOfElements = size(T,1);
nOfElementNodes = size(T,2);
coordRef_star = referenceElement_star.NodesCoord;
npoints = size(coordRef_star,1);
nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);

% Vandermonde matrix
coordRef = referenceElement.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = referenceElement.degree;
V = Vandermonde_LP(nDeg,coordRef);
%invV = inv(V');

% % Compute shape functions at interpolation points
% shapeFunctions = zeros(npoints,nOfNodes);
% for ipoint = 1:npoints
%     p = orthopoly2D(coordRef_star(ipoint,:),nDeg);
%     shapeFunctions(ipoint,:) = (V'\p)';
% end

%Shape functions for interpolation of nodal values
[N,dNdxi,dNdeta]=evaluateNodalBasisTri(coordRef_star,coordRef,nDeg);


% u star initialization
u_star = zeros(2*npoints*nOfElements,1);


% Loop in elements
for ielem = 1:nOfElements
    
    c=find(ielem==Elements.D1);
    d=find(ielem==Elements.D2);
    
    Te_lin=T(ielem,:);
    Xold = X(Te_lin,:);
    Xg=N*Xold;
    LSe_star = EvaluateLS(Xg,1);
    
    ind = (ielem-1)*(2*nOfElementNodes) + (1:2*nOfElementNodes);
    u1 = N*u(ind(1:nOfElementNodes));
    u2 = N*u(ind(nOfElementNodes+1:end));
    
    ue=[u1;u2];
    
    x_ind=(2*ind-1);
    y_ind=(2*ind);
    
    qe(1:2:2*npoints,1) = N*q(x_ind(1:nOfElementNodes));
    qe(2:2:2*npoints,1) = N*q(y_ind(1:nOfElementNodes));
    
    qe(2*npoints+1:2:4*npoints-1,1) =  N*q(x_ind((nOfElementNodes+1:end)));
    qe(2*npoints+2:2:4*npoints,1)   =  N*q(y_ind(nOfElementNodes+1:end));
    
    
    [Ke,Beqe,int_ue_star, int_ue,int_u_star_H,int_u_H] = ElementalMatrices(Xg,referenceElement_star,ue,qe,Elements,ielem,LSe_star,mu1,mu2);      
     
    if isempty(c) && isempty(d);
        
        K1 = [Ke int_ue_star; int_ue_star' 0];
        f1 = [Beqe;int_ue];
    
        K2 = [K1 [int_u_star_H;0] ; int_u_star_H' 0 0];
        f2 = [f1;int_u_H];
        sol=K2\f2;
        
        sol=[sol(1:end-2)];
        
    else
        
    K = [Ke((1:nOfElementNodes_star),(1:nOfElementNodes_star)) int_ue_star(1:nOfElementNodes_star); int_ue_star(1:nOfElementNodes_star)' 0];
    f = [Beqe(1:nOfElementNodes_star); int_ue];
        
    sol=K\f;
    sol=[sol(1:(end-1)) ; zeros(nOfElementNodes_star,1)];        
        
    end
    
    index = (ielem-1)*(2*npoints) + (1:2*npoints);
    % postprocessed solution
    u_star(index) = sol;
end

%%
%% ELEMENTAL MATRIX

function [K,Bq,int_u_star,int_u,int_u_star_H,int_u_H] = ElementalMatrices(Xe,referenceElement_star,ue,qe,Elements,iElem,LSe_star,mu1,mu2)


c=find(iElem==Elements.D1);
d=find(iElem==Elements.D2);
nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);

size_matrix=2*nOfElementNodes_star;


K = zeros(size_matrix,size_matrix);
Bq = zeros(size_matrix,1);
int_u_star = zeros(size_matrix,1);
int_u = 0 ;
int_u_star_H=0;
int_u_H =0;


if ~isempty(c)
    
    % Information of the reference element
    IPw = referenceElement_star.IPweights;
    Nold = referenceElement_star.N;
    Nxi_old = referenceElement_star.Nxi;
    Neta_old = referenceElement_star.Neta;
    N=[Nold, Nold*0];
    Nxi=[Nxi_old,Nxi_old*0];
    Neta=[Neta_old,Neta_old*0];
    ngauss=size(IPw,1);
    weight=IPw;
    mu=ones(ngauss,1)*mu1;  % !!!!!! ATTENTION HERE !!!
    H=ones(ngauss,1);
    
elseif ~isempty(d)
    
    % Information of the reference element
    IPw = referenceElement_star.IPweights;
    Nold = referenceElement_star.N;
    Nxi_old = referenceElement_star.Nxi;
    Neta_old = referenceElement_star.Neta;
    N=[Nold, Nold*0];
    Nxi=[Nxi_old,Nxi_old*0];
    Neta=[Neta_old,Neta_old*0];
    ngauss=size(IPw,1);
    weight=IPw;
    mu=ones(ngauss,1)*mu2;  % !!!!!! ATTENTION HERE !!!
    H=ones(ngauss,1);
    
else isempty(c) && isempty(d);
    
    
    p_star = referenceElement_star.degree;
    %Quadrature for standart triangle and quadrilateral
    [zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*p_star);
    wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
    [zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*p_star);
    
    
    [zgp,wgp,n1,n2,PtsInt,FaceInfo] = ModifyQuadrature(LSe_star,referenceElement_star,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement_star.degree,referenceElement_star.NodesCoord,zgp);
    Nold = shapeFunctions(:,:,1)';  Nxi_old = shapeFunctions(:,:,2)'; Neta_old= shapeFunctions(:,:,3)';
    weight=wgp;
    ngauss=n1+n2;
    mu=[(mu1.*ones(n1,1));(mu2*ones(n2,1))];
    H=[(1.*ones(n1,1));(-1*ones(n2,1))];
    
    %     %Enrichment of Shape Functions
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
    
    
end


% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% VOLUME COMPUTATIONS: LOOP IN GAUSS POINTS
for g = 1:ngauss
    
    %Shape functions and derivatives at the current integration point
    Nold_g=Nold(g,:);
    Nxi_g_old = Nxi_old(g,:);
    Neta_g_old = Neta_old(g,:);
    
    N_g=N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    % gauss pts coordinates
    xg=Nold_g*[xe ye];
    
    %Jacobian
    J = [Nxi_g_old*xe	  Nxi_g_old*ye
        Neta_g_old*xe     Neta_g_old*ye];
    
    %Integration weight
    dvolu=weight(g)*det(J);
    
    %x and y derivatives
    invJ = inv(J);
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    
    % u and q at gauss points
    u_g=N_g*ue;
    qx_g=N_g*qe(1:2:end);
    qy_g =N_g*qe(2:2:end);
    
    
    %Contribution of the current integration point to the elemental matrix
    Bq = Bq - (Nx_g'*qx_g + Ny_g'*qy_g)*dvolu;
    K = K + (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu*mu(g);
    int_u_star = int_u_star + N_g'*dvolu;
    int_u = int_u + u_g*dvolu;
    int_u_star_H = int_u_star_H + H(g)*N_g'*dvolu;
    int_u_H = int_u_H + H(g)*u_g*dvolu;
end

