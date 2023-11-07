
function [KK,f, QQ, UU, Qf, Uf] = hdg_matrix(X,T,F,referenceElement,infoFaces,tau,mu,u0,q0)
% This routine is valid only for triangles

nOfElements = size(T,1);
nOfElementNodes = size(T,2);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(infoFaces.intFaces,1)+size(infoFaces.extFaces,1);
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
    indu0=(iElem-1)*(nOfElementNodes)+(1:nOfElementNodes);
    u_0=u0(indu0)';
    indq0=(iElem-1)*(2*nOfElementNodes)+(1:2*nOfElementNodes);
    q_0=q0(indq0)';

    %u_0=Xe(:,1)';
    %q_0=[-1  0   -1  0   -1  0];
    %uh_0= [ 0  0  0  0   0   0];

    % elemental matrices
    [Qe,Ue,Qfe,Ufe,Alq,Alu,All,Rl] = KKeElementalMatrices(Xe,referenceElement,tau,iElem,mu,u_0,q_0);
    
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
    ffe = Rl-(Alq*Qfe + Alu*Ufe);
    
    %All faces are assembled (Dirichlet included)
    aux = (1:nOfFaceNodes);
    ind = [(Fe(1)-1)*nOfFaceNodes + aux,(Fe(2)-1)*nOfFaceNodes + aux,(Fe(3)-1)*nOfFaceNodes + aux];
    KK(ind,ind)=KK(ind,ind)+ KKe;
    f(ind) = f(ind) + ffe;
end



% ELEMENTAL MATRIX
function [Q,U,Qf,Uf,Alq,Alu,All,Rl] = KKeElementalMatrices(Xe,referenceElement,tau,iElem,mu,u_0,q_0)

nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)
%tau = [tau,zeros(1,nOfFaces-1)]; %Non-null only for 1st face in each element

fe = zeros(nOfElementNodes,1);
Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
f_Aqqx = zeros(nOfElementNodes,1);
f_Aqqy = zeros(nOfElementNodes,1);
f_Aqq = zeros(2*nOfElementNodes,1);
Auq = zeros(nOfElementNodes,2*nOfElementNodes);
f_Auq = zeros(nOfElementNodes,1);
Aqu = zeros(2*nOfElementNodes,nOfElementNodes);
f_Aqu = zeros(2*nOfElementNodes,1);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
Aql = zeros(2*nOfElementNodes,3*nOfFaceNodes);
f_Aql = zeros(2*nOfElementNodes,1);
f_Alq = zeros(3*nOfFaceNodes,1);
Auu = zeros(nOfElementNodes,nOfElementNodes);
f_Auu_s = zeros(nOfElementNodes,1);
f_Auu_v = zeros(nOfElementNodes,1);
f_Auu_2 = zeros(nOfElementNodes,1);
Alu = zeros(3*nOfFaceNodes,nOfElementNodes);
Aul = zeros(nOfElementNodes,3*nOfFaceNodes);
All = zeros(3*nOfFaceNodes,3*nOfFaceNodes);
NN_g = zeros(2,2*nOfElementNodes); %Shape functions for vector variable q
NN_xy = zeros(1,2*nOfElementNodes);
f_All = zeros(3*nOfFaceNodes,1);
f_Alu_2 = zeros(3*nOfFaceNodes,1);
Rw = zeros(2*nOfElementNodes,1);
Fv = zeros(nOfElementNodes,1);
Rl = zeros(3*nOfFaceNodes,1);

% Information of the reference element
N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;  
N = referenceElement.N;
Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
IPw = referenceElement.IPweights;
ngauss=size(IPw,1);
weight=IPw;

%Numerical quadrature
Xg = N*Xe; %physical coordinates
u0_g = N*u_0';
q0x_g = N*q_0(1:2:end)';
q0y_g = N*q_0(2:2:end)';
q0x=q_0(1:2:end);
q0y=q_0(2:2:end);

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
    dvolu=weight(g)*det(J);
    
    %x and y derivatives
    invJ = inv(J);
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    NN_xy(1:2:end) = Nx_g;    NN_xy(2:2:end) = Ny_g;
    
    %%% Contribution of the current integration point to the elemental matrix
    Aqq = Aqq + NN_g'*NN_g*dvolu/mu;  % Aqq
    Auq = Auq + N_g'*NN_xy*dvolu;     % Auq
    Aqu=  Aqu - NN_xy'*N_g*dvolu;
    grad_v = [Nx_g ; Ny_g];
    Auu = Auu - ([u0_g(g)*2, u0_g(g)*2]*grad_v)'*N_g*dvolu; %Auu
    %Auu = Auu -  ([u0_g(g); u0_g(g)].*N_g)'*grad_v*dvolu; %Auu
   
    %source term.
    fe = fe + N_g'*sourceTerm(g)*dvolu;


    %%%%% RHS terms coming from linearization %%%%%%%%
    f_Aqqx = f_Aqqx - q0x_g(g).*N_g'*dvolu/mu;
    f_Aqqy = f_Aqqy - q0y_g(g).*N_g'*dvolu/mu;
    f_Aqu = f_Aqu + u0_g(g).*NN_xy'*dvolu;
    f_Auq = f_Auq - (q0x*Nx_g' + q0y*Ny_g').*N_g'*dvolu;
    f_Auu_v = f_Auu_v + ([u0_g(g)*u0_g(g)  u0_g(g)*u0_g(g)]*grad_v)'*dvolu;
end
    f_Aqq(1:2:end,1) = f_Aqqx;
    f_Aqq(2:2:end,1) = f_Aqqy;

%% Faces computations:
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces

    %tau_f = tau(iface);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    xf = xe(nodes);
    yf = ye(nodes);
    uh_0face=u_0(nodes);
    %uh_0face=uh_0(nodes);
    q_0fx=q_0(1:2:end);
    q_0fx= q_0fx(nodes);
    q_0fy=q_0(2:2:end);
    q_0fy= q_0fy(nodes);
   
    % Inizialization
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    Alq_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    Alq_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    Aql_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    Aql_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    f_Aql_fnx = zeros(nOfFaceNodes,1);
    f_Aql_fny = zeros(nOfFaceNodes,1);
    f_Alq_fn = zeros(nOfFaceNodes,1);
    %f_Alq_fny = zeros(nOfFaceNodes,1);
    Auu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    f_Auu_f = zeros(nOfFaceNodes,1);
    f_Auu_f2 = zeros(nOfFaceNodes,1);
    f_Auu_f3 = zeros(nOfFaceNodes,1);
    Aul_f = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_f2 = zeros(nOfFaceNodes,nOfFaceNodes);
    All_f = zeros(nOfFaceNodes);
    All_f2 = zeros(nOfFaceNodes);
    
    
    glim=ngauss_f;
    %Shape functions and derivatives
    Nf_g = N1d;
    Nfxi_g = Nx1d;
    %Integration weight
    IPwFace=IPw_f;
    uh_0faceg=Nf_g*uh_0face';
    q_0fxg=Nf_g*q_0fx';
    q_0fyg=Nf_g*q_0fy';

    %  LOOP IN GAUSS POINTS
    for g = 1:glim
        
        % Integration weight
        xyDer_g = Nfxi_g(g,:)*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline=IPwFace(g)*xyDerNorm_g;
        
        % Unit normal to the boundary
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) -t_g(1)];
             
        tau_f=1;            
            %%%Contribution of the current integration point to the elemental matrix
            Aql_fnx = Aql_fnx + (Nf_g(g,:)*n_g(1))'*Nf_g(g,:)*dline;
            Aql_fny = Aql_fny + (Nf_g(g,:)*n_g(2))'*Nf_g(g,:)*dline;            
            Alq_fnx = Alq_fnx + Nf_g(g,:)'*Nf_g(g,:)*n_g(1)*dline;
            Alq_fny = Alq_fny + Nf_g(g,:)'*Nf_g(g,:)*n_g(2)*dline;
            Auu_f = Auu_f + tau_f*Nf_g(g,:)'*Nf_g(g,:)*dline; 
            Aul_f = Aul_f - tau_f*Nf_g(g,:)'*Nf_g(g,:)*dline; 
            %Aul_f2 = Aul_f2 + [uh_0faceg(g)  uh_0faceg(g)]*n_g'*Nf_g(g,:)'*Nf_g(g,:)*dline;
            %Aul_f2 = Aul_f2 + [uh_0faceg(g) ; uh_0faceg(g)].*Nf_g(g,:)*n_g'*Nf_g(g,:)*dline;
            All_f = All_f -tau_f*Nf_g(g,:)'*Nf_g(g,:)*dline;
            All_f2 = All_f2 + [uh_0faceg(g)*2  uh_0faceg(g)*2]*n_g'*Nf_g(g,:)'*Nf_g(g,:)*dline;

            %%%%%%%%%%%%%%%  RHS vectors  %%%%%%%%%%%%%%%%%%%% 
            f_Aql_fnx = f_Aql_fnx - (Nf_g(g,:)*n_g(1))'*uh_0faceg(g)*dline;
            f_Aql_fny = f_Aql_fny - (Nf_g(g,:)*n_g(2))'*uh_0faceg(g)*dline;
            f_Auu_f = f_Auu_f - ([uh_0faceg(g)*uh_0faceg(g) ;  uh_0faceg(g)*uh_0faceg(g)]'*n_g'*Nf_g(g,:))'*dline;
            f_Auu_f2 = f_Auu_f2 - tau_f*uh_0faceg(g).*Nf_g(g,:)'*dline; 
            f_Auu_f3 = f_Auu_f3 + tau_f*uh_0faceg(g).*Nf_g(g,:)'*dline; 
            f_Alq_fn = f_Alq_fn - Nf_g(g,:)'*[q_0fxg(g), q_0fyg(g)]*n_g'*dline;
            %f_Alq_fny = f_Alq_fny - Nf_g(g,:)'*q_0fyg(g)*n_g(2)*dline;
   
    end
        
    Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Alq_fnx;
    Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Alq_fny;
    Aql(2*nodes-1,ind_face) = Aql(2*nodes-1,ind_face) + Aql_fnx;
    Aql(2*nodes,ind_face) = Aql(2*nodes,ind_face) + Aql_fny;
    Auu(nodes,nodes) = Auu(nodes,nodes) + Auu_f; 
    Aul(nodes,ind_face) = Aul(nodes,ind_face) + Aul_f +All_f2; 
    Alu(ind_face,nodes) = Alu(ind_face,nodes) + Auu_f;
    All(ind_face,ind_face) = All(ind_face,ind_face) + All_f + All_f2;

    %%%%%%%%%%%%%%%%%%%%%  RHS vectors %%%%%%%%%%%%%%%%%%%%%%
 
    f_Aql(2*nodes-1,1) = f_Aql(2*nodes-1,1) + f_Aql_fnx;
    f_Aql(2*nodes,1) = f_Aql(2*nodes,1) + f_Aql_fny;
    f_Auu_s(nodes,1) = f_Auu_s(nodes,1) + f_Auu_f; 
    f_Auu_2(nodes,1) = f_Auu_2(nodes,1) + f_Auu_f2 + f_Auu_f3;
    f_Alq(ind_face,1) = f_Alq(ind_face,1) + f_Alq_fn;
    f_All (ind_face,1)= f_All(ind_face,1) + f_Auu_f; 
    f_Alu_2(ind_face,1) = f_Alu_2(ind_face,1) + f_Auu_f2 + f_Auu_f3;

end

    Rw = Rw+f_Aqq+f_Aqu+f_Aql; 
    Fv = Fv+fe+f_Auq+f_Auu_s+f_Auu_v+f_Auu_2; 
    Rl = Rl+f_Alq+f_All+f_Alu_2;


    % Elemental mapping
    A = [Auu Auq; Aqu Aqq];
    UQ = -A\[Aul;Aql];
    fUQ= A\[Fv;Rw];
    U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
    Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q


     %%% Test
     u_analy=analiticalSolutionLaplace(Xe);
     x = Xe(:,1);
     y = Xe(:,2);  
     q_analy=zeros(2*length(u_analy),1);
     du_dx = 0+2.*x;
     du_dy = 0+0.*y;
     q_analy(1:2:end)=-mu.*du_dx;
     q_analy(2:2:end)=-mu.*du_dy;

     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama_x=Xe(ce,1);
     coord_gama_y=Xe(ce,2);
     coord_gama=[coord_gama_x,coord_gama_y];
     gama=analiticalSolutionLaplace(coord_gama);    

     test_result_1=Fv-((Auu*u_analy)+(Auq*q_analy)+(Aul*gama));
     test_result_2=Rw-((Aqq*q_analy)+(Aqu*u_analy)+(Aql*gama));
     %test_result_3=Rl-((Alq*q_analy)+(Alu*u_analy)+(All*gama));

     disp(['Standard Element: Test result equation 1 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_1)));
     disp(['Standard Element: Test result equation 2 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_2)));














