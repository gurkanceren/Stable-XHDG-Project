function [KK,f, QQ, UU, Qf, Uf] = HDGMatrixStraightSides(muElem,X,T,F,referenceElement,infoFaces,tau)
% [KK,f, QQ, UU, Qf, Uf] = HDGMatrixStraightSides(muElem,X,T,F,referenceElement,infoFaces,tau)
% This routine is valid only for triangles
% Elemental computations HDG for Laplace, for straight-sided elements

%Computation of elemental matrices in the reference element
[Me,Cxi,Ceta,Mf] = ElementalMatricesLaplaceReferenceElement(referenceElement);

nOfFaces = max(max(F));
nOfElements = size(T,1);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nDOF = nOfFaces*nOfFaceNodes;
%KK = spalloc(nDOF,nDOF,nOfFaceNodes*nDOF);
f = zeros(nDOF,1);
QQ = cell(nOfElements,1); UU = cell(nOfElements,1); Qf = cell(nOfElements,1); Uf = cell(nOfElements,1);

% loop in elements
indK = 1; indf=1; n = nOfElements*(3*nOfFaceNodes)^2; ind_i  = zeros(1,n); ind_j  = zeros(1,n); coef_K = zeros(1,n);
for iElem = 1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)

    % elemental matrices
    [Qe,Ue,Qfe,Ufe,Alq,Alu,All] = KKeElementalMatricesStraightSides(muElem(iElem),Xe,referenceElement,tau(iElem,:),Me,Cxi,Ceta,Mf);
    
    % Interior faces seen from the second element are flipped to have
    % proper orientation
    flipFace = zeros(1,3); %Boolean (1=face to be flipped)
    flipFace(isFeInterior)=any(infoFaces.intFaces(Fe(isFeInterior),[1 3])<iElem,2);
    indL=1:3*nOfFaceNodes; 
    aux=nOfFaceNodes:-1:1; indflip=[aux,nOfFaceNodes+aux,2*nOfFaceNodes+aux];
    aux=ones(1,nOfFaceNodes); aux = [flipFace(1)*aux, flipFace(2)*aux, flipFace(3)*aux];
    indL(aux==1)=indflip(aux==1); %permutation for local numbering
    
    Qe=Qe(:,indL);    Ue=Ue(:,indL);
    Alq=Alq(indL,:);  Alu=Alu(indL,:);   All=All(indL,indL);
    
    %The local problem solver is stored for postprocess
    QQ{iElem} = sparse(Qe);  UU{iElem} = sparse(Ue);
    Qf{iElem} = sparse(Qfe); Uf{iElem} = sparse(Ufe);
    
    %Elemental matrices to be assembled
    KKe = Alq*Qe + Alu*Ue + All;
    ffe = -(Alq*Qfe + Alu*Ufe);
    
    %All faces are assembled (Dirichlet included)
    aux = (1:nOfFaceNodes);
    indRC = [(Fe(1)-1)*nOfFaceNodes + aux,(Fe(2)-1)*nOfFaceNodes + aux,(Fe(3)-1)*nOfFaceNodes + aux];
    f(indRC) = f(indRC) + ffe;
    %KK(indRC,indRC)=KK(indRC,indRC)+ KKe;
    for irow = 1:(3*nOfFaceNodes)
        for icol = 1:(3*nOfFaceNodes)
            ind_i(indK)  = indRC(irow); ind_j(indK)  = indRC(icol);
            coef_K(indK) = KKe(irow,icol); indK = indK+1;
        end
    end
end
KK = sparse(ind_i,ind_j,coef_K);


%%
%% ELEMENTAL MATRICES straight-sided element
function [Q,U,Qf,Uf,Alq,Alu,All] = KKeElementalMatricesStraightSides(mu,Xe,referenceElement,tau,Me,Cxi,Ceta,Mf)

nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

%% Volume integrals
%Jacobian for straight-sided triangle
J = [(Xe(2,:)-Xe(1,:))/2 ; (Xe(3,:)-Xe(1,:))/2];
invJ= inv(J);
detJ = det(J);
Me = detJ*Me;

%Element matrices
Aqq = zeros(2*size(Me));
aux = 1:2:2*nOfElementNodes; aux2 = 2:2:2*nOfElementNodes;
Aqq(aux,aux)=Me; Aqq(aux2,aux2)=Me;
Auq = zeros(nOfElementNodes,2*nOfElementNodes);
Auq(:,aux) =detJ*( invJ(1,1)*Cxi + invJ(1,2)*Ceta ); %x derivatives & 1st component of q
Auq(:,aux2)=detJ*( invJ(2,1)*Cxi + invJ(2,2)*Ceta ); %y derivatives & 2nd component of q

%Faces matrices
Auu = zeros(nOfElementNodes,nOfElementNodes);
Alu = zeros(3*nOfFaceNodes,nOfElementNodes);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
All = zeros(3*nOfFaceNodes,3*nOfFaceNodes);
for iface = 1:nOfFaces
    tau_f = tau(iface);
    % Nodes in the face
    nodes = faceNodes(iface,:);
    Xf = Xe(nodes,:);
    %Jacobian of face transformation
    vec_t = Xf(end,:)-Xf(1,:);
    dline = norm(vec_t)/2;
    vec_t = vec_t/norm(vec_t);
    vec_n = [vec_t(2) -vec_t(1)];
    %Face matrices
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    Alq(ind_face,2*nodes-1) = Mf*(vec_n(1)*dline);
    Alq(ind_face,2*nodes) = Mf*(vec_n(2)*dline);
    Auu_f = Mf*(tau_f*dline);  
    Auu(nodes,nodes) = Auu(nodes,nodes) + Auu_f;
    Alu(ind_face,nodes) = Alu(ind_face,nodes) + Auu_f;
    All(ind_face,ind_face) = -Auu_f;
end

%Computation of r.h.s. source term (analytical laplacian)
Xg = N*Xe;
sourceTerm = -mu*analiticalLaplacianLaplace(Xg);
fe = N'*(detJ*sourceTerm.*(referenceElement.IPweights));


% Elemental mapping
Aqu = -mu*Auq'; Aul = -Alu'; Aql = mu*Alq';
A = [Auu Auq; Aqu Aqq];
UQ = -A\[Aul;Aql];
fUQ= A\[fe;zeros(2*nOfElementNodes,1)];
U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q


save matricesStraightSides A Aul Aql All fe sourceTerm

%%________________________________________________________________________
%% ELEMENTAL MATRICES for reference element
function [Me,Cxi,Ceta,Mf] = ElementalMatricesLaplaceReferenceElement(referenceElement)

% Information of the reference element
N = referenceElement.N;
Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
N1d = referenceElement.N1d; 
%Numerical quadrature
IPw = referenceElement.IPweights; ngauss=size(IPw,1);
IPw_f = referenceElement.IPweights1d;

%% Volume computation: "loop" in Gauss points
ngauss = length(IPw);
dvolu = spdiags(IPw,0,ngauss,ngauss);
Me = N'*(dvolu*N);
Cxi = N'*(dvolu*Nxi);
Ceta = N'*(dvolu*Neta);
%% Face coputations
ngauss_f = length(IPw_f);
Mf = N1d'*(spdiags(IPw_f',0,ngauss_f,ngauss_f)*N1d);



%% OLD VERSION
% % nOfElementNodes = size(referenceElement.NodesCoord,1);
% % nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
% % 
% % %Inicialization of matrices
% % Me = zeros(nOfElementNodes,nOfElementNodes);
% % Cxi = zeros(nOfElementNodes,nOfElementNodes);
% % Ceta = zeros(nOfElementNodes,nOfElementNodes);
% % Mf = zeros(nOfFaceNodes,nOfFaceNodes);
% % 
%% Volume computation: loop in Gauss points
% % for g = 1:ngauss
% %     %Shape functions and derivatives at the current integration point
% %     N_g = N(g,:);
% %     Nxi_g = Nxi(g,:);    Neta_g = Neta(g,:);
% %     %Integration weight
% %     dvolu=IPw(g);
% %     %Contribution of the current integration point to the elemental matrix
% %     Me = Me + N_g'*N_g*dvolu; %mass matrix
% %     Cxi = Cxi + N_g'*Nxi_g*dvolu;% N'*dN_dxi
% %     Ceta = Ceta + N_g'*Neta_g*dvolu;% N'*dN_deta
% % end
% % %% Face coputations
% % ngauss_f = length(IPw_f);
% % %  LOOP IN GAUSS POINTS
% % for g = 1:ngauss_f
% %     Nf_g = N1d(g,:);
% %     Mf = Mf + Nf_g'*Nf_g*IPw_f(g);
% % end


