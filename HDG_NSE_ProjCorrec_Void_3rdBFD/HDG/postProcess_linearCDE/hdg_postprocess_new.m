function [u_star,shapeFunctions] = hdg_postprocess_new(X,T,u,q,qt,referenceElement_star,referenceElement,mu_vector)

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

% qt star initialization
qt_star = zeros(2*npoints*nOfElements,1);

% Loop in elements
for iElem = 1:nOfElements
    %qt_star = zeros(2*npoints,1);
    
    Te_lin = T(iElem,:);
    Xold = X(Te_lin,:);
    Xe=shapeFunctions*Xold;
    mu=mu_vector(iElem);

    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    ind_star = (iElem-1)*npoints+1:iElem*npoints;
    ind_starq = (iElem-1)*2*npoints+1:iElem*2*npoints;
    %ind_star2 = (iElem-1)*2*npoints+1:iElem*2*npoints;

    % projection of u & u_hat into the new polynomial space.
    ue = shapeFunctions*u(ind);
    %ue_new = 0.*ue;
    %ue_hat = shapeFunctions*uhat(ind);

    % projection of qe & qte into the new polynomial space.
    qe(1:2:2*npoints,1) = shapeFunctions*q(2*ind-1);
    qe(2:2:2*npoints,1) = shapeFunctions*q(2*ind);
   
    qte(1:2:2*npoints,1) = shapeFunctions*qt(2*ind-1);
    qte(2:2:2*npoints,1) = shapeFunctions*qt(2*ind);

    %qte_new = 0.*qte;

    %% postprocessed solution for qt

    %[K1,Bq_1,K2,Bq_2] = ElementalMatricesforq(Xe,referenceElement_star,qte);

    %Kq = [K1 ;K2];
    %fq = [Bq_1;Bq_2];

    %sol_q = Kq\fq;

    %qt_star(ind_starq) = sol_q;

    %qe = qe';
    %q_x = qe(1:2:2*npoints);
    %q_y = qe(2:2:2*npoints);
    %qte = qte';
    %qt_x = qte(1:2:2*npoints);
    %qt_y = qe(2:2:2*npoints);

    %% postprocessed solution for u

    [Ke1,Beqe1,Ke2,Beqe2,Ke3,Beqe3,int_ue_star,int_ue] = ElementalMatricesforu(Xe,referenceElement_star,ue,qte,mu);
    
    % Lagrange multipliers

    %K = [Ke int_ue_star; int_ue_star' 0];
    %f = [Beqe;int_ue];
    %rr = length(int_ue_star);
    %rs = zeros(rr,1);

    %K = [Ke2 int_ue_star; int_ue_star' 0]; 
    %f = [Beqe2;int_ue];
    

    K = [Ke1 int_ue_star; Ke2 int_ue_star; Ke3 int_ue_star; int_ue_star' 0]; 
    f = [Beqe1;Beqe2;Beqe3;int_ue];

    sol = K\f;
        
    u_star(ind_star) = sol(1:end-1);

    %disp('Hola');
end
%disp('Hola');

%%
%% ELEMENTAL MATRIX for u

function [K_1,Bq1,K_2,Bq2,K_3,Bq3,int_u_star,int_u] = ElementalMatricesforu(Xe,referenceElement_star,ue,qte,mu)

%**************************************************************
nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);
nOfFaceNodes = size(referenceElement_star.NodesCoord1d,1);
nOfFaces = size(referenceElement_star.faceNodes,1);
faceNodes = referenceElement_star.faceNodes;
%**************************************************************
% Initialization
Bq1 = zeros(nOfElementNodes_star,1);
Bq2 = zeros(nOfElementNodes_star,1);
Bq3 = zeros(nOfElementNodes_star,1);
Bq_surf = zeros(nOfElementNodes_star,1);
int_u_star = zeros(nOfElementNodes_star,1);
int_u = 0;

%K = zeros(nOfElementNodes_star,nOfElementNodes_star);
Kvol1 = zeros(nOfElementNodes_star,nOfElementNodes_star);
Kvol2 = zeros(nOfElementNodes_star,nOfElementNodes_star);
Kvol3 = zeros(nOfElementNodes_star,nOfElementNodes_star);
Ku = zeros(nOfElementNodes_star,nOfElementNodes_star);
Kustar = zeros(nOfElementNodes_star,nOfElementNodes_star);
Ksurf = zeros(nOfElementNodes_star,nOfElementNodes_star);
Kq = zeros(nOfElementNodes_star,nOfElementNodes_star);
Kcu = zeros(nOfElementNodes_star,nOfElementNodes_star);
%**************************************************************
% Information of the reference element
IPw = referenceElement_star.IPweights;
N = referenceElement_star.N;
Nxi = referenceElement_star.Nxi;
Neta = referenceElement_star.Neta;
NN_xy = zeros(1,2*nOfElementNodes_star);
NN_g = zeros(2,2*nOfElementNodes_star); %Shape functions for vector variable q
IPw_f = referenceElement_star.IPweights1d;
N1d = referenceElement_star.N1d; 
Nx1d = referenceElement_star.N1dxi;
%**************************************************************

% Number of Gauss points in the interior
ngauss = length(IPw);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

l_d = 1; % representative diffusive length scale
c_x = 1; % advection speed of the signal in x-direction
c_y = 1; % advection speed of the signal in y-direction

%pp_q = length(qe);
%qe_x = qe(1:2:pp_q); 
%qe_y = qe(2:2:pp_q);

pp_qt = length(qte);
qte_x = qte(1:2:pp_qt);
qte_y = qte(2:2:pp_qt);

% call the subroutine which computes the source term.
sourceTerm = @analiticalLaplacianLaplace;

%**************************************************************
%% VOLUME COMPUTATIONS: LOOP IN GAUSS POINTS
%**************************************************************
for g = 1:ngauss
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    NN_g(1,1:2:end)=N_g; NN_g(2,2:2:end)=N_g;
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    
    %Integration weight
    dvolu=IPw(g)*det(J);
    
    %x and y derivatives
    invJ = inv(J); 
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    NN_xy(1:2:end) = Nx_g;
    NN_xy(2:2:end) = Ny_g;
    
    % u and grad(qt) at gauss points
    u_g = N_g*ue;
    qtx_g = Nx_g*qte_x; 
    qty_g = Ny_g*qte_y; 
    
    qt_gx = N_g*qte_x;
    qt_gy = N_g*qte_y;
    
    %Contribution of the current integration point to the elemental matrix
    ba = (N_g'*qtx_g + N_g'*qty_g)*dvolu;
    Bq2 = Bq2 + ba;

    bb = (Nx_g'*qt_gx + Ny_g'*qt_gy)*dvolu;
    Bq3 = Bq3 + (ba+bb);

    a = (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu/mu;
    b = (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    Kvol1 = Kvol1 + (a+b);
    %bs = (N_g'*Nx_g + N_g'*Ny_g)*dvolu;
    %Kvol1 = Kvol1 + (a-bs);

    c = (N_g'*Nx_g + N_g'*Ny_g)*dvolu;
    d = (N_g'*Nx_g + N_g'*Ny_g)*dvolu;
    Kvol2 = Kvol2 + (c+d);
    %cs = (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    %ds = (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    %Kvol2 = Kvol2 - (cs+ds);
    
    e = (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    f = (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    g = (N_g'*Nx_g + N_g'*Ny_g)*dvolu;
    h = (N_g'*Nx_g + N_g'*Ny_g)*dvolu;
    Kvol3 = Kvol3 + (e+f+g+h);

    int_u_star = int_u_star + N_g'*dvolu;
    int_u = int_u + u_g*dvolu;
end

K_1 = Kvol1; K_2 = Kvol2; K_3 = Kvol3;

%disp('Hola');
%{
%
%**************************************************************
%% FACES COMPUTATIONS:
%**************************************************************
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces
    
    % tau_f = tau(iface);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    xf = xe(nodes);
    yf = ye(nodes);

    u_f = ue(nodes);
    %uhat_f = ue_hat(nodes);
    %q_fx = qe_x(nodes);
    %q_fy = qe_y(nodes);
    qt_fx = qte_x(nodes);
    qt_fy = qte_y(nodes);
    
    % Gauss points position
    xfg = N1d*xf;
    yfg = N1d*yf;   
  
    % Inizialization
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    Kq_fx = zeros(nOfFaceNodes,nOfFaceNodes);
    Kq_fy = zeros(nOfFaceNodes,nOfFaceNodes);
    Kcu_fx = zeros(nOfFaceNodes,nOfFaceNodes);
    Kcu_fy = zeros(nOfFaceNodes,nOfFaceNodes);
    Ku_f = zeros(nOfFaceNodes,nOfFaceNodes);
    Ku_starf = zeros(nOfFaceNodes,nOfFaceNodes);
    Bq_fg = zeros(nOfFaceNodes,1);
    
    %  LOOP IN GAUSS POINTS
    for g = 1:ngauss_f
        
        % Shape functions and derivatives at the current integration point
        Nf_g = N1d(g,:);
        Nfxi_g = Nx1d(g,:);
        
        % Integration weight
        xyDer_g = Nfxi_g*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline=IPw_f(g)*xyDerNorm_g;
        
        % Unit normal to the boundary
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) -t_g(1)];

        % compute convective-tau and diffusive-tau.
        tau_cf = abs(c_x.*n_g(1)+c_y.*n_g(2));  % |c.n|
        tau_df = (mu/l_d);  %*tau_f;
        % compute combined tau.
        tau_f = tau_cf+tau_df;
        
        %
        % u and q at the gauss points;
        u_gf = Nf_g*u_f;
        %uhat_gf = Nf_g*uhat_f;
        %q_gfx = Nf_g*q_fx;
        %q_gfy = Nf_g*q_fy;
        qt_gfx = Nf_g*qt_fx;
        qt_gfy = Nf_g*qt_fy;

        % Contribution of the current integration point to the elemental matrix
        Bq_fg = Bq_fg + Nf_g'*(qt_gfx*n_g(1)+qt_gfy*n_g(2))*dline;
        
        Kq_fx = Kq_fx + Nf_g'*Nf_g*n_g(1)*dline;
        Kq_fy = Kq_fy + Nf_g'*Nf_g*n_g(2)*dline;

        Kcu_fx = Kcu_fx + Nf_g'*Nf_g*n_g(1)*dline;
        Kcu_fy = Kcu_fy + Nf_g'*Nf_g*n_g(2)*dline;

        Ku_f = Ku_f + tau_f*(Nf_g')*Nf_g*dline;
        Ku_starf = Ku_starf - tau_f*(Nf_g')*Nf_g*dline;
 
    end
    Bq_surf(nodes) = Bq_surf(nodes) + Bq_fg;

    Ku(nodes,nodes) = Ku(nodes,nodes) + Ku_f;
    Kustar(nodes,nodes) = Kustar(nodes,nodes) + Ku_starf;

    Kq(nodes,nodes) = Kq(nodes,nodes) + (Kq_fx+Kq_fy);
    Kcu(nodes,nodes) = Kcu(nodes,nodes) + (Kcu_fx+Kcu_fy);

end

Ksurf = Kq+Kcu+Ku+Kustar;

K_1 = Kvol1 + Kcu;
K_2 = Kvol2 + Ksurf;
K_3 = Ksurf;

Bq3 = Bq_surf;
%}

%%
%% ELEMENTAL MATRIX for qt

function [K1,Bq_1,K2,Bq_2] = ElementalMatricesforq(Xe,referenceElement_star,qte)

%**************************************************************
nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);
nOfFaceNodes = size(referenceElement_star.NodesCoord1d,1);
nOfFaces = size(referenceElement_star.faceNodes,1);
faceNodes = referenceElement_star.faceNodes;
%**************************************************************
% Initialization
Bqx_1 = zeros(nOfElementNodes_star,1);
Bqy_1 = zeros(nOfElementNodes_star,1);
Bqx_2 = zeros(nOfElementNodes_star,1);
Bqy_2 = zeros(nOfElementNodes_star,1);
Bq_1 = zeros(2*nOfElementNodes_star,1);
Bq_2 = zeros(2*nOfElementNodes_star,1);
%Bq3 = zeros(nOfElementNodes_star,1);
%Bq_surf = zeros(nOfElementNodes_star,1);
%int_u_star = zeros(nOfElementNodes_star,1);
%int_u = 0;

%K = zeros(nOfElementNodes_star,nOfElementNodes_star);
Kvol1 = zeros(2*nOfElementNodes_star,2*nOfElementNodes_star);
Kvol2 = zeros(2*nOfElementNodes_star,2*nOfElementNodes_star);
%Kvol3 = zeros(nOfElementNodes_star,nOfElementNodes_star);
%Ku = zeros(nOfElementNodes_star,nOfElementNodes_star);
%Kustar = zeros(nOfElementNodes_star,nOfElementNodes_star);
%Ksurf = zeros(nOfElementNodes_star,nOfElementNodes_star);
%Kq = zeros(nOfElementNodes_star,nOfElementNodes_star);
%Kcu = zeros(nOfElementNodes_star,nOfElementNodes_star);
%**************************************************************
% Information of the reference element
IPw = referenceElement_star.IPweights;
N = referenceElement_star.N;
Nxi = referenceElement_star.Nxi;
Neta = referenceElement_star.Neta;
NN_xy = zeros(1,2*nOfElementNodes_star);
NN_g = zeros(1,2*nOfElementNodes_star); %Shape functions for vector variable q
IPw_f = referenceElement_star.IPweights1d;
N1d = referenceElement_star.N1d; 
Nx1d = referenceElement_star.N1dxi;
%**************************************************************

% Number of Gauss points in the interior
ngauss = length(IPw);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

l_d = 1; % representative diffusive length scale
c_x = 1; % advection speed of the signal in x-direction
c_y = 1; % advection speed of the signal in y-direction

%pp_q = length(qe);
%qe_x = qe(1:2:pp_q); 
%qe_y = qe(2:2:pp_q);

pp_qt = length(qte);
qte_x = qte(1:2:pp_qt);
qte_y = qte(2:2:pp_qt);

%**************************************************************
%% VOLUME COMPUTATIONS: LOOP IN GAUSS POINTS
%**************************************************************
for g = 1:ngauss
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    NN_g(1:2:end)=N_g; NN_g(2:2:end)=N_g;
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    
    %Integration weight
    dvolu=IPw(g)*det(J);
    
    %x and y derivatives
    invJ = inv(J); 
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    NN_xy(1:2:end) = Nx_g;
    NN_xy(2:2:end) = Ny_g;
    
    % u and grad(qt) at gauss points
    %u_g = N_g*ue;
    qtx_g = Nx_g*qte_x; 
    qty_g = Ny_g*qte_y; 
    
    qt_gx = N_g*qte_x;
    qt_gy = N_g*qte_y;
    
    %Contribution of the current integration point to the elemental matrix
    ba_x = (N_g'*qtx_g)*dvolu;
    ba_y = (N_g'*qty_g)*dvolu;

    bb_x = (Nx_g'*qt_gx)*dvolu;
    bb_y = (Ny_g'*qt_gy)*dvolu;

    Bqx_1 = Bqx_1 + (ba_x+bb_x);
    Bqy_1 = Bqy_1 + (ba_x+bb_x);

    Bq_1(1:2:end) = Bqx_1;
    Bq_1(2:2:end) = Bqy_1;

    a = (NN_xy'*NN_xy)*dvolu;
    b = (NN_g'*NN_xy)*dvolu;
    Kvol1 = Kvol1 + (a+b);

    Bqx_2 = Bqx_2 + bb_x;
    Bqy_2 = Bqy_2 + bb_y;

    Bq_2(1:2:end) = Bqx_2;
    Bq_2(2:2:end) = Bqy_2;

    Kvol2 = Kvol2 + (NN_xy'*NN_xy)*dvolu;

    %e = (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    %f = (Nx_g'*Nx_g + Ny_g'*Ny_g)*dvolu;
    %g = (N_g'*Nx_g + N_g'*Ny_g)*dvolu;
    %h = (N_g'*Nx_g + N_g'*Ny_g)*dvolu;
    %Kvol3 = Kvol3 + (e+f+g+h);

    %int_u_star = int_u_star + N_g'*dvolu;
    %int_u = int_u + u_g*dvolu;
end

K1 = Kvol1; K2 = Kvol2; %K_3 = Kvol3;

%disp('Hola');