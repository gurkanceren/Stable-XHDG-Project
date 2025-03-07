function [u_star,shapeFunctions] = hdg_postprocess(X,T,u,q,referenceElement_star,referenceElement,Elements,mu1,mu2)

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
u_star = zeros(2*npoints*nOfElements,1);


% Loop in elements
for ielem = 1:nOfElements
    
    
    c=find(ielem==Elements.D1);
    d=find(ielem==Elements.D2);
    
    Te_lin=T(ielem,:);
    Xold = X(Te_lin,:);
    Xg=shapeFunctions*Xold;
    LSe_star = EvaluateLS(Xg,2);
    
    if isempty(c) && isempty(d);  %cut element
        
        ind = (ielem-1)*(2*nOfElementNodes) + (1:2*nOfElementNodes);
        u1 = shapeFunctions*u(ind(1:nOfElementNodes));
        u2 = shapeFunctions*u(ind(nOfElementNodes+1:end));
        
        ue=[u1;u2];
        
        x_ind=(2*ind-1);
        y_ind=(2*ind);
        
        qe(1:2:2*npoints,1) = shapeFunctions*q(x_ind(1:nOfElementNodes));
        qe(2:2:2*npoints,1) = shapeFunctions*q(y_ind(1:nOfElementNodes));
        
        qe(2*npoints+1:2:4*npoints-1,1) = shapeFunctions *q(x_ind((nOfElementNodes+1:end)));
        qe(2*npoints+2:2:4*npoints,1)   =  shapeFunctions*q(y_ind(nOfElementNodes+1:end));
        
    else  %standard element
        
        ind = (ielem-1)*(2*nOfElementNodes) + (1:2*nOfElementNodes);
        ue = shapeFunctions*u(ind(1:nOfElementNodes));
        qe(1:2:2*npoints,1) = shapeFunctions*q(2*ind(1:nOfElementNodes)-1);
        qe(2:2:2*npoints,1) = shapeFunctions*q(2*ind(1:nOfElementNodes));
    end
    
    
    
    [Ke,Beqe,int_ue_star, int_ue] = ElementalMatrices(Xg,referenceElement_star,ue,qe,Elements,ielem,LSe_star,referenceElement,mu1,mu2);
    
    % Lagrange multipliers
    K = [Ke int_ue_star; int_ue_star' 0];
    f = [Beqe;int_ue];
    
    % elemental solution
    sol = K\f;
    
    if isempty(c) && isempty(d);
        
         sol=sol(1:end-1);
         index = (ielem-1)*(2*npoints) + (1:2*npoints);
        
    else
        
        sol=[sol(1:end-1);zeros(npoints,1)];
        index = (ielem-1)*(2*npoints) + (1:2*npoints);
    end
    
    % postprocessed solution
    u_star(index) = sol;
    
    qe=[];
end

%%
%% ELEMENTAL MATRIX

function [K,Bq,int_u_star,int_u] = ElementalMatrices(Xe,referenceElement_star,ue,qe,Elements,iElem,LSe_star,referenceElement,mu1,mu2)


c=find(iElem==Elements.D1);
d=find(iElem==Elements.D2);

nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);

if ~isempty(c) || ~isempty(d)
    
    size_matrix=nOfElementNodes_star;
    
else isempty(c) && isempty(d);
    
    size_matrix=2*nOfElementNodes_star;
end




K = zeros(size_matrix,size_matrix);
Bq = zeros(size_matrix,1);
int_u_star = zeros(size_matrix,1);
int_u = 0;


if ~isempty(c)
    
    % Information of the reference element
    IPw = referenceElement_star.IPweights;
    Nold = referenceElement_star.N;
    Nxi_old = referenceElement_star.Nxi;
    Neta_old = referenceElement_star.Neta;
    
    N=Nold;
    Nxi=Nxi_old;
    Neta=Neta_old ;
    ngauss=size(IPw,1);
    weight=IPw;
    mu=mu1;  % !!!!!! ATTENTION HERE !!!
    
elseif ~isempty(d)
    
    % Information of the reference element
    IPw = referenceElement_star.IPweights;
    Nold = referenceElement_star.N;
    Nxi_old = referenceElement_star.Nxi;
    Neta_old = referenceElement_star.Neta;
    
    N=Nold;
    Nxi=Nxi_old;
    Neta=Neta_old ;
    ngauss=size(IPw,1);
    weight=IPw;
    mu=mu2;  % !!!!!! ATTENTION HERE !!!
    
    
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
    
    
end


% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% VOLUME COMPUTATIONS: LOOP IN GAUSS POINTS
for g = 1:ngauss
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    Nold_g=Nold(g,:);
    Nxi_g_old = Nxi_old(g,:);
    Neta_g_old = Neta_old(g,:);
    
    
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
    %     NN_xy(1:2:end) = Nx_g;
    %     NN_xy(2:2:end) = Ny_g;
    
    % u and q at gauss points
    u_g =  N_g*ue;
    qx_g = N_g*qe(1:2:end);
    qy_g = N_g*qe(2:2:end);
    
    
    % viscosity (for cut elements)
    if  g <= (ngauss/2) && isempty(c) && isempty(d)       
        mu=mu1;
    elseif g > (ngauss/2) && isempty(c) && isempty(d)
        mu=mu2;       
    else
        mu=mu;
    end

    
    %Contribution of the current integration point to the elemental matrix
    Bq = Bq - (Nx_g'*qx_g + Ny_g'*qy_g)*dvolu;
    K = K + (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu*mu;
    int_u_star = int_u_star + N_g'*dvolu;
    int_u = int_u + u_g*dvolu;
end


