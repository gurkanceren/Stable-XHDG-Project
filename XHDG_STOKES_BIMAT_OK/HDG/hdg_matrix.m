function [KK,f, QQ, UU, Qf, Uf] = hdg_matrix(X,T,F,referenceElement,infoFaces,tau)
% This routine is valid only for triangles

nOfElements = size(T,1);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = max(max(F));
nDOF = nOfFaces*nOfFaceNodes;
KK = spalloc(nDOF,nDOF,nOfFaceNodes*nDOF);
f = zeros(nDOF,1);
QQ = cell(nOfElements,1); UU = cell(nOfElements,1); Qf = cell(nOfElements,1); Uf = cell(nOfElements,1);

% loop in elements
for iElem = 1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)

    % elemental matrices
    [Qe,Ue,Qfe,Ufe,Alq,Alu,All] = KKeElementalMatricesLaplace(Xe,referenceElement,tau);
    
    % Interior faces seen from the second element are flipped to have
    % proper orientation
    flipFace = zeros(1,3); %Boolean (1=face to be flipped)
    flipFace(isFeInterior)=any(infoFaces.intFaces(Fe(isFeInterior),[1 3])<iElem,2);
    indL=1:3*nOfFaceNodes; 
    aux=nOfFaceNodes:-1:1; indflip=[aux,nOfFaceNodes+aux,2*nOfFaceNodes+aux];
    aux=ones(1,nOfFaceNodes); aux = [flipFace(1)*aux, flipFace(2)*aux, flipFace(3)*aux];
    indL(aux==1)=indflip(aux==1); %permutation for local numbering
    
    Qe=Qe(:,indL);
    Ue=Ue(:,indL);
    Alq=Alq(indL,:);
    Alu=Alu(indL,:);
    All=All(indL,indL);
    
    %The local problem solver is stored for postprocess
    QQ{iElem} = sparse(Qe);  UU{iElem} = sparse(Ue);
    Qf{iElem} = sparse(Qfe); Uf{iElem} = sparse(Ufe);
    
    %Elemental matrices to be assembled
    KKe = Alq*Qe + Alu*Ue + All;
    ffe = -(Alq*Qfe + Alu*Ufe);
    
    %All faces are assembled (Dirichlet included)
    aux = (1:nOfFaceNodes);
    ind = [(Fe(1)-1)*nOfFaceNodes + aux,(Fe(2)-1)*nOfFaceNodes + aux,(Fe(3)-1)*nOfFaceNodes + aux];
    KK(ind,ind)=KK(ind,ind)+ KKe;
    f(ind) = f(ind) + ffe;
end



%%
%% ELEMENTAL MATRIX
function [Q,U,Qf,Uf,Alq,Alu,All] = KKeElementalMatricesLaplace(Xe,referenceElement,tau)

nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)
%tau = [tau,zeros(1,nOfFaces-1)]; %Non-null only for 1st face in each element

fe = zeros(nOfElementNodes,1);
Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
Auq = zeros(nOfElementNodes,2*nOfElementNodes);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
Auu = zeros(nOfElementNodes,nOfElementNodes);
Alu = zeros(3*nOfFaceNodes,nOfElementNodes);
All = zeros(3*nOfFaceNodes,3*nOfFaceNodes);
NN_g = zeros(2,2*nOfElementNodes); %Shape functions for vector variable q
NN_xy = zeros(1,2*nOfElementNodes);

% Information of the reference element
N = referenceElement.N;
Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;
%Numerical quadrature
IPw = referenceElement.IPweights; ngauss=size(IPw,1);
Xg = N*Xe; %physical coordinates
IPw_f = referenceElement.IPweights1d;

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% analytical laplacian
sourceTerm = analiticalLaplacianLaplace(Xg);

%% Volume computation: loop in Gauss points
for g = 1:ngauss
    
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    NN_g(1,1:2:end)=N_g; NN_g(2,2:2:end)=N_g;
    Nxi_g = Nxi(g,:);    Neta_g = Neta(g,:);
    
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end

    %Integration weight
    dvolu=IPw(g)*det(J);
    
    %x and y derivatives
    invJ = inv(J);  
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    NN_xy(1:2:end) = Nx_g;    NN_xy(2:2:end) = Ny_g;
    
    %Contribution of the current integration point to the elemental matrix
    Aqq = Aqq + NN_g'*NN_g*dvolu;
    Auq = Auq + N_g'*NN_xy*dvolu;
    fe = fe - N_g'*sourceTerm(g)*dvolu;
end
inv_Aqq = inv(Aqq);

%% Faces coputations:
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces
    
    tau_f = tau(iface);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    xf = xe(nodes);
    yf = ye(nodes);
    
    % Gauss points position
    xfg = N1d*xf;
    yfg = N1d*yf;   
  
    % Inizialization
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    Alq_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    Alq_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    Auu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    All_f = zeros(nOfFaceNodes);
    
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
                
        % Contribution of the current integration point to the elemental matrix
        Alq_fnx = Alq_fnx + Nf_g'*Nf_g*n_g(1)*dline;
        Alq_fny = Alq_fny + Nf_g'*Nf_g*n_g(2)*dline;
        Auu_f = Auu_f + tau_f *Nf_g'*Nf_g*dline;
        All_f = All_f - tau_f *Nf_g'*Nf_g*dline;
    end
    Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Alq_fnx;
    Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Alq_fny;
    Auu(nodes,nodes) = Auu(nodes,nodes) + Auu_f;
    Alu(ind_face,nodes) = Alu(ind_face,nodes) + Auu_f;
    All(ind_face,ind_face) = All(ind_face,ind_face) + All_f;
end

% Elemental mapping
Aqu = -Auq'; Aul = -Alu'; Aql = Alq';
A = [Auu Auq; Aqu Aqq];
UQ = -A\[Aul;Aql];
fUQ= A\[fe;zeros(2*nOfElementNodes,1)];
U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q

save matrices A Aul Aql All fe sourceTerm

