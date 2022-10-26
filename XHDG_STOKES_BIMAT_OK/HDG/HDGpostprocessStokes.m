function [u_star,N] = HDGpostprocessStokes(X,T,u,L,referenceElement_star,referenceElement,Elements,example)

nOfElements = size(T,1);
coordRef=referenceElement.NodesCoord;
nOfElementNodes = size(coordRef,1);
coordRef_star = referenceElement_star.NodesCoord;
npoints = size(coordRef_star,1);
degree=referenceElement.degree;
nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);

%Shape functions for interpolation of nodal values
[N,dNdxi,dNdeta]=evaluateNodalBasisTri(coordRef_star,coordRef,degree);

% u star initialization
u_star = zeros(2*npoints*nOfElements,2);

% Loop in elements
for iElem = 1:nOfElements   
    ind = (iElem-1)*2*nOfElementNodes+1:iElem*2*nOfElementNodes;
    ind_star = (iElem-1)*2*npoints+1:iElem*2*npoints;
    
    Te = T(iElem,:);
    Xe = X(Te ,:);
    Xe_star=N*Xe;
    LSe_star = EvaluateLS(Xe_star,example);   %p+1 degree LS function
    %1st component
    %ue = N*u(ind,1);
    
    u1 = N*u(ind(1:nOfElementNodes),1);
    u1H = N*u(ind(nOfElementNodes+1:end),1); 
    ue=[u1;u1H];
    
    qe(1:2:2*npoints,1) = N*L(ind(1:nOfElementNodes),1);
    qe(2:2:2*npoints,1) = N*L(ind(1:nOfElementNodes),2);
    qe(2*npoints+1:2:4*npoints-1,1) = N*L(ind(nOfElementNodes+1:end),1);
    qe(2*npoints+2:2:4*npoints,1)   = N*L(ind(nOfElementNodes+1:end),2);
    
    c=find(iElem==Elements.D1); d=find(iElem==Elements.D2);
        
    [Ke,Beqe,int_ue_star, int_ue,int_u_star_H,int_u_H] = ElementalMatrices(X(T(iElem,:),:),Xe_star,referenceElement_star,ue,qe,Elements,LSe_star,iElem);
    if isempty(c) && isempty(d);
        
        K1 = [Ke int_ue_star; int_ue_star' 0];
        f1 = [Beqe;int_ue];
        
        K2 = [K1 [int_u_star_H;0] ; int_u_star_H' 0 0];
        f2 = [f1;int_u_H];
        sol=K2\f2;
        u_star(ind_star,1) =sol(1:end-2);
        
    else
        
    K = [Ke((1:nOfElementNodes_star),(1:nOfElementNodes_star)) int_ue_star(1:nOfElementNodes_star); int_ue_star(1:nOfElementNodes_star)' 0];
    f = [Beqe(1:nOfElementNodes_star); int_ue];
        
    sol=K\f;
    sol=[sol(1:(end-1)) ; zeros(nOfElementNodes_star,1)];        
    u_star(ind_star,1) = sol;    
    end
    

    
    %2nd component
    u2 = N*u(ind(1:nOfElementNodes),2);
    u2H = N*u(ind(nOfElementNodes+1:end),2); 
    ue=[u2;u2H];
    qe(1:2:2*npoints,1) = N*L(ind(1:nOfElementNodes),3);
    qe(2:2:2*npoints,1) = N*L(ind(1:nOfElementNodes),4);
    qe(2*npoints+1:2:4*npoints-1,1) = N*L(ind(nOfElementNodes+1:end),3);
    qe(2*npoints+2:2:4*npoints,1)   = N*L(ind(nOfElementNodes+1:end),4);
    
    [Ke,Beqe,int_ue_star, int_ue,int_u_star_H,int_u_H] = ElementalMatrices(X(T(iElem,:),:),Xe_star,referenceElement_star,ue,qe,Elements,LSe_star,iElem);
    if isempty(c) && isempty(d);
        
        K1 = [Ke int_ue_star; int_ue_star' 0];
        f1 = [Beqe;int_ue];
    
        K2 = [K1 [int_u_star_H;0] ; int_u_star_H' 0 0];
        f2 = [f1;int_u_H];
        sol=K2\f2;
        u_star(ind_star,2) =sol(1:end-2);
        
    else
        
    K = [Ke((1:nOfElementNodes_star),(1:nOfElementNodes_star)) int_ue_star(1:nOfElementNodes_star); int_ue_star(1:nOfElementNodes_star)' 0];
    f = [Beqe(1:nOfElementNodes_star); int_ue];
        
    sol=K\f;
    sol=[sol(1:(end-1)) ; zeros(nOfElementNodes_star,1)];        
    u_star(ind_star,2) = sol;    
    end

    
end

%%
%% ELEMENTAL MATRIX

function [K,Bq,int_u_star,int_u,int_u_star_H,int_u_H] = ElementalMatrices(Xe,Xe_star,referenceElement_star,ue,qe,Elements,LSe_star,iElem)
global mu1 mu2;
nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);
K = zeros(2*nOfElementNodes_star,2*nOfElementNodes_star);
Bq = zeros(2*nOfElementNodes_star,1);
int_u_star = zeros(2*nOfElementNodes_star,1);
int_u = 0;
int_u_star_H=0;
int_u_H=0;

c=find(iElem==Elements.D1); d=find(iElem==Elements.D2);

if ~isempty(c)
    
% Information of the reference element
IPw = referenceElement_star.IPweights;
Nold = referenceElement_star.N;
Nxi_old = referenceElement_star.Nxi;
Neta_old = referenceElement_star.Neta;
N=[Nold, Nold*0];
Nxi=[Nxi_old,Nxi_old*0];
Neta=[Neta_old,Neta_old*0];
NN_xy = zeros(1,2*nOfElementNodes_star);
% Number of Gauss points in the interior
ngaussend = length(IPw);
weight = IPw;
H=ones(ngaussend,1);

elseif ~isempty(d)
    
% Information of the reference element
IPw = referenceElement_star.IPweights;
Nold = referenceElement_star.N;
Nxi_old = referenceElement_star.Nxi;
Neta_old = referenceElement_star.Neta;
N=[Nold, Nold*0];
Nxi=[Nxi_old,Nxi_old*0];
Neta=[Neta_old,Neta_old*0];
NN_xy = zeros(1,2*nOfElementNodes_star);
% Number of Gauss points in the interior
ngaussend = length(IPw);
weight = IPw;
H=ones(ngaussend,1);

else  isempty(c) && isempty(d);
    
    p_star = referenceElement_star.degree;
    %Quadrature for standart triangle and quadrilateral
    [zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*p_star);
    wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
    [zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*p_star);
    
    
    [zgp,wgp,n1,n2,PtsInt,FaceInfo] = ModifyQuadrature(LSe_star,referenceElement_star,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement_star.degree,referenceElement_star.NodesCoord,zgp);
    Nold = shapeFunctions(:,:,1)';  Nxi_old = shapeFunctions(:,:,2)'; Neta_old= shapeFunctions(:,:,3)';
    weight=wgp; ngaussend=n1+n2;
    H=[(1.*ones(n1,1));(-1*ones(n2,1))];
    
    Nl=Nold(1:n1,:)*1;
    Nr=Nold(n1+1:n2+n1,:)*-1;
    NH=[Nl;Nr];
    N=[Nold,NH];
    
    Nlxi=Nxi_old(1:n1,:)*1;
    Nrxi=Nxi_old(n1+1:n2+n1,:)*-1;
    NHxi=[Nlxi;Nrxi];
    Nxi=[Nxi_old,NHxi];
    
    Nleta=Neta_old(1:n1,:)*1;
    Nreta=Neta_old(n1+1:n2+n1,:)*-1;
    NHeta=[Nleta;Nreta];
    Neta=[Neta_old,NHeta];
    
end

% x and y coordinates of the element nodes
xe = Xe_star(:,1); ye = Xe_star(:,2);

% VOLUME COMPUTATIONS: LOOP IN GAUSS POINTS
for g = 1:ngaussend
    
    Nold_g=Nold(g,:);
    Nxi_g_old = Nxi_old(g,:);
    Neta_g_old = Neta_old(g,:);
    
    N_g=N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    xy_g=Nold_g*[xe ye];
    
    %Jacobian
    J = [Nxi_g_old*xe	  Nxi_g_old*ye
        Neta_g_old*xe     Neta_g_old*ye];
    
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
    qx_g = N_g*qe(1:2:end);
    qy_g = N_g*qe(2:2:end);
    
    %Contribution of the current integration point to the elemental matrix
    Bq = Bq + (Nx_g'*qx_g + Ny_g'*qy_g)*dvolu; 
    K = K + (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    int_u_star = int_u_star + N_g'*dvolu;
    int_u = int_u + u_g*dvolu;
    int_u_star_H = int_u_star_H + H(g)*N_g'*dvolu;
    int_u_H = int_u_H + H(g)*u_g*dvolu;
end


