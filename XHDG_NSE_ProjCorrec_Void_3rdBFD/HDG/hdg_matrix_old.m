
function [KK,f, QQ, UU, Qf, Uf] = hdg_matrix(X,T,F,referenceElement,infoFaces,tau,LS,Elements,mu,c_x,c_y,upwind,dt,cnt,u,time)
% This routine is valid only for triangles

nOfElements = size(T,1);
nOfElementNodes = size(referenceElement.NodesCoord,1);
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
    LSe = LS(Te);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)
    u_Elem = u((iElem-1)*nOfElementNodes+(1:nOfElementNodes));
    
    % elemental matrices
    [Qe,Ue,Qfe,Ufe,Alq,Alu,All] = KKeElementalMatrices(Xe,referenceElement,tau,LSe,iElem,Elements,mu,c_x,c_y,upwind,dt,cnt,u,time,u_Elem);
    
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



% ELEMENTAL MATRIX
function [Q,U,Qf,Uf,Alq,Alu,All] = KKeElementalMatrices(Xe,referenceElement,tau,LSe,iElem,Elements,mu,c_x,c_y,upwind,dt,cnt,u,time,u_Elem)


c=find(iElem==Elements.D1);
d=find(iElem==Elements.D2);

%c_x = 1;
%c_y = 1;
l_d = 1;

nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)
%tau = [tau,zeros(1,nOfFaces-1)]; %Non-null only for 1st face in each element

fe_s = zeros(nOfElementNodes,1);
fe_uns = zeros(nOfElementNodes,1);
Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
Auq = zeros(nOfElementNodes,2*nOfElementNodes);
%Aqu = zeros(2*nOfElementNodes,nOfElementNodes);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
Auu = zeros(nOfElementNodes,nOfElementNodes);
Auu_vol_s = zeros(nOfElementNodes,nOfElementNodes); % Auu_vol = (cu)*grad(w)
Auu_vol_t = zeros(nOfElementNodes,nOfElementNodes);
Auu_surf = zeros(nOfElementNodes,nOfElementNodes);
Alu = zeros(3*nOfFaceNodes,nOfElementNodes);
Aul = zeros(nOfElementNodes,3*nOfFaceNodes);
All = zeros(3*nOfFaceNodes,3*nOfFaceNodes);
NN_g = zeros(2,2*nOfElementNodes); %Shape functions for vector variable q
NN_xy = zeros(1,2*nOfElementNodes);

% Information of the reference element
N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;



if ~isempty(c)
    
    N = referenceElement.N;
    Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
    IPw = referenceElement.IPweights;
    ngauss=size(IPw,1);
    weight=IPw;
    PtsInt=[];
    FaceInfo=[];
    
    
elseif ~isempty(d)
    
    N=0;
    ngauss=0;
    PtsInt=[];
    FaceInfo=[];
    
else isempty(c) && isempty(d);
    
    
    p = referenceElement.degree;
    %Quadrature for standart triangle and quadrilateral
    [zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*p);
    wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
    [zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*p);
    
    
    [zgp,wgp,n1,n2,PtsInt,FaceInfo] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zgp);
    N = shapeFunctions(:,:,1)';  Nxi = shapeFunctions(:,:,2)'; Neta= shapeFunctions(:,:,3)';
    weight=wgp; ngauss=n1;   %SFM

    
    
end


%Numerical quadrature
Xg = N*Xe; %physical coordinates
IPw_f = referenceElement.IPweights1d;

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% analytical laplacian
sourceTerm = analiticalLaplacianLaplace(Xg,c_x,c_y,time);

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
    xy_g=N_g*Xe;
    
    %x and y derivatives
    invJ = inv(J);
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    NN_xy(1:2:end) = Nx_g;    NN_xy(2:2:end) = Ny_g;
    
    %%% Contribution of the current integration point to the elemental matrix
    Aqq = Aqq + NN_g'*NN_g*dvolu/mu;  % Aqq
    Auq = Auq + N_g'*NN_xy*dvolu;     % Auq
   % Aqu=  Aqu + NN_g(1,:)'*N_g*dvolu;

   % computing of Auu_vol.
    c_vec = [c_x c_y];
    cu = c_vec'*N_g; %cu_t = cu';
    grad_w = [Nx_g ; Ny_g];
    grad_w_t = grad_w';
    Auu_vol_temp = -grad_w_t*cu;
    Auu_vol_s = Auu_vol_s + Auu_vol_temp*dvolu; %
    Auu_vol_t = Auu_vol_t + (N_g'*N_g*dvolu/dt);

   % vol. integral of the source term.
    %fe = fe - N_g'*sourceTerm(g)*dvolu;
    fe_s = fe_s + N_g'*sourceTerm(g)*dvolu;

    if(cnt==1)
        u_old_gauss = initial_condition(xy_g);
        %u_old_gauss = analiticalSolutionLaplace(xy_g,time);
    else
        u_old_gauss = N_g*u_Elem; 
    end 
    %u_old_gauss = N_g*u_Elem;
   
    fe_uns = fe_uns + N_g'*u_old_gauss*dvolu/dt;

end


%% Faces computations:
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces
    
    
    if ~isempty(FaceInfo)
        k=find(FaceInfo(:,1)==iface);
    else
        k=[];
    end
    
    %tau_f = tau(iface);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    xf = xe(nodes);
    yf = ye(nodes);
    
    % Inizialization
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    Alq_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    Alq_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    Auu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_f = zeros(nOfFaceNodes,nOfFaceNodes);
    All_f = zeros(nOfFaceNodes);
    
    
    
    if isempty(k) && any(LSe(faceNodes(iface,:))>0);  % NORMAL FACES
        
        glim=ngauss_f;
        %Shape functions and derivatives
        Nf_g = N1d;
        Nfxi_g = Nx1d;
        %Integration weight
        IPwFace=IPw_f;
        
    elseif ~isempty(k);    %CUT FACE
        
        % Shape Functions for cut faces:
        
        [zgp_f,wgp_f,n1_f,n2_f,IntPt] = ModifyQuadrature1D(LSe(faceNodes(iface,:),1),referenceElement);
        shapeFunctions_f=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord1d,zgp_f);
        glim=n1_f;
        
        % Shape functions and derivatives
        Nf_g =shapeFunctions_f(:,:,1)';
        Nfxi_g = shapeFunctions_f(:,:,2)';
        %Integration weight
        IPwFace=wgp_f;
        
        
    else %VOID FACE
        
        glim=0;
        
    end
    
    
    %  LOOP IN GAUSS POINTS
    for g = 1:glim
        
        % Integration weight
        xyDer_g = Nfxi_g(g,:)*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline=IPwFace(g)*xyDerNorm_g;
        
        % Unit normal to the boundary
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) -t_g(1)];
        
        % compute convective and diffusive eta:
        eta_c = abs(c_x.*n_g(1)+c_y.*n_g(2));  % |c.n|
        eta_d = (mu/l_d);

        if (upwind==1)
            % compute c.n+ and |c.n+|:
            a = abs(c_x.*n_g(1)+c_y.*n_g(2));  % |c.n+|
            b = c_x.*n_g(1)+c_y.*n_g(2);       %  c.n+
            % compute tau_convective:
            if (b==0)
                tau_cf = eta_c;
            else
                tau_cf = eta_c*(a+b)/(2*a);
            end
            % compute tau_diffusive:
            if (b==0)
                tau_df = eta_d;
            else
                tau_df = eta_d*(a+b)/(2*a);
            end 
            %
        else
            %
            % compute convective-tau and diffusive-tau.
            tau_cf = eta_c; %abs(c_x.*n_g(1)+c_y.*n_g(2));  % |c.n|
            tau_df = eta_d; %(mu/l_d);  
            %
        end

        % compute combined tau.
        tau_f = tau_cf+tau_df;
        % compute (c.n) as pp.
        pp = c_x.*n_g(1)+c_y.*n_g(2);
        % compute (c.n-tau) as bb.
        if (pp>0)
            bb = -tau_df;
        elseif (pp==0)
            bb = -tau_df;
        else
            bb = pp-tau_f;
        end

        %bb = 0-tau_f;
        aa = 0.5*pp; %aa = (1/2)*(c.n)
         
        if(tau_f>=aa)
            %
            % Contribution of the current integration point to the elemental matrix
            Alq_fnx = Alq_fnx + Nf_g(g,:)'*Nf_g(g,:)*n_g(1)*dline;
            Alq_fny = Alq_fny + Nf_g(g,:)'*Nf_g(g,:)*n_g(2)*dline;
            Auu_f = Auu_f + tau_f*Nf_g(g,:)'*Nf_g(g,:)*dline; % Auu_surf
            Aul_f = Aul_f + bb*Nf_g(g,:)'*Nf_g(g,:)*dline; 
            All_f = All_f + bb*Nf_g(g,:)'*Nf_g(g,:)*dline;
        else
            disp('Lemma 3.1 violated');
        end    

    end
    
    
    Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Alq_fnx;
    Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Alq_fny;
    Auu_surf(nodes,nodes) = Auu_surf(nodes,nodes) + Auu_f; %
    Aul(nodes,ind_face) = Aul(nodes,ind_face) + Aul_f; %
    Alu(ind_face,nodes) = Alu(ind_face,nodes) + Auu_f;
    All(ind_face,ind_face) = All(ind_face,ind_face) + All_f;
    
end


if  isempty(c) && isempty (d)
    
    %Calculation of new AquuhatI,AuuhatI and AuuI matrices because of applied
    %Dirichlet boundary conditions on I
    
    p=referenceElement.degree;
    g=length(referenceElement.IPweights1d);
    AquhatI=zeros(2*nOfElementNodes,nOfFaceNodes);
    
    
    PtsInt=reshape(PtsInt',2,p+1)'; %interface p+1 nodes in the REFERENCE element
    
    zg=referenceElement.N1d*PtsInt;%integration points on the interface in the REFERENCE element
    
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zg);
    Nzg = shapeFunctions(:,:,1)';  % 2d shape functions on the refence element at integration points
    
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,PtsInt);
    NPtsInt = shapeFunctions(:,:,1)';  %2D Shape functions on the REFERENCE element at interface nodes
    NPtsInt_1 = NPtsInt';
    
    Pphy=NPtsInt*Xe; %p+1 interface nodes on the REAL element
    Pgauss = referenceElement.N1d*Pphy;
       
    Iprime=referenceElement.N1dxi*Pphy;   %*Pphy;
    % I=referenceElement.N1d*Pphy;
    %dxdxi = Nx1d*xf; dydxi = Nx1d*yf;
    dxdxiNorm = sqrt(Iprime(:,1).^2+Iprime(:,2).^2);
    
    %Initialization
    
    AuuhatI=zeros(nOfElementNodes,nOfFaceNodes);
    AuuI=zeros(nOfElementNodes,nOfElementNodes);
    AquhatI_x_f=zeros(nOfElementNodes,nOfFaceNodes);
    AquhatI_y_f=zeros(nOfElementNodes,nOfFaceNodes);
    
    for igauss=1:g
        
        
        % Integration weight
        normIprime = norm(Iprime(igauss,:));
        wg=referenceElement.IPweights1d(igauss)*normIprime;
        
        % Unit normal to the boundary
        t_g = Iprime(igauss,:)/normIprime;
        n_g = [t_g(2) -t_g(1)];
        
        % compute convective and diffusive eta:
        eta_c = abs(c_x.*n_g(1)+c_y.*n_g(2));  % |c.n|
        eta_d = (mu/l_d);

        if (upwind==1)
            % compute c.n+ and |c.n+|:
            a = abs(c_x.*n_g(1)+c_y.*n_g(2));  % |c.n+|
            b = c_x.*n_g(1)+c_y.*n_g(2);       %  c.n+
            % compute tau_convective:
            if (b==0)
                tau_cf = eta_c;
            else
                tau_cf = eta_c*(a+b)/(2*a);
            end
            % compute tau_diffusive:
            if (b==0)
                tau_df = eta_d;
            else
                tau_df = eta_d*(a+b)/(2*a);
            end 
            %
        else
            %
            % compute convective-tau and diffusive-tau.
            tau_cf = eta_c; %abs(c_x.*n_g(1)+c_y.*n_g(2));  % |c.n|
            tau_df = eta_d; %(mu/l_d);  
            %
        end

        % compute combined tau.
        tau_f = tau_cf+tau_df;
        % compute (c.n) as pp.
        pp = c_x.*n_g(1)+c_y.*n_g(2);
        % compute (c.n-tau) as bb.
        if (pp>0)
            bb = -tau_df;
        elseif (pp==0)
            bb = -tau_df;
        else
            bb = pp-tau_f;
        end

        %{
        AuuhatI = AuuhatI - tau_f*Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*wg;
        AuuI = AuuI + tau_f*Nzg(igauss,:)'*Nzg(igauss,:)*wg;
        AquhatI_x_f =  AquhatI_x_f + Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(1)*wg;
        AquhatI_y_f =  AquhatI_y_f + Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(2)*wg;
        %}

        AuuhatI = AuuhatI - bb*Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*wg;
        AuuI = AuuI + tau_f*Nzg(igauss,:)'*Nzg(igauss,:)*wg;
        AquhatI_x_f =  AquhatI_x_f - Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(1)*wg;
        AquhatI_y_f =  AquhatI_y_f - Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(2)*wg;

    end
    
    AquhatI(1:2:end,:) = AquhatI_x_f;
    AquhatI(2:2:end,:) = AquhatI_y_f;
    
    
    uhatI_e_old=analiticalSolutionLaplace(Pphy,time);
    uhatI_e=analiticalSolutionLaplace(Pgauss,time);
    dline = spdiags(dxdxiNorm.*IPw_f,0,g,g);
    M = (referenceElement.N1d)'*(dline*referenceElement.N1d);
    b = (referenceElement.N1d)'*(dline*uhatI_e);    
    uface = M\b;

    
    Auu=Auu_vol_t+Auu_vol_s+Auu_surf+AuuI;
    
    %     figure(1); hold on;
    %     plot(Pphy(:,1),Pphy(:,2),'ro');
    %     quiver(I(:,1),I(:,2),Iprime(:,1),Iprime(:,2));
    
    fu_I = (AuuhatI*uface);
    fq_I = (AquhatI*uface);
    
    %Elemental mapping
    
    Aqu = -Auq'; 
    %Aul = -Alu'; 
    Aql = Alq';
    A = [Auu Auq; Aqu Aqq];
    UQ = -A\[Aul;Aql];
    %fUQ= A\[(fe-(AuuhatI*uhatI_e));-(AquhatI*uhatI_e)];
    fUQ= A\[(fe_uns+fe_s+fu_I);fq_I];
    U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
    Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q
        
   %{
    % %     %%%TEST  TEST
     u_analy=analiticalSolutionLaplace(Xe,time);
     x = Xe(:,1);
     y = Xe(:,2);
%     
     q_analy=zeros(2*length(u_analy),1);
     du_dx = 0+y.*y;
     du_dy = 0+2.*x.*y;
     q_analy(1:2:end)=-1.*du_dx;
     q_analy(2:2:end)=-1.*du_dy;
%     
%         
     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama_x=Xe(ce,1);
     coord_gama_y=Xe(ce,2);
     coord_gama=[coord_gama_x,coord_gama_y];
     gama=analiticalSolutionLaplace(coord_gama,time);
%     
     uhat_analy= analiticalSolutionLaplace(Pphy,time);
     a1 = Aqu*u_analy;
     a2 = Aqq*q_analy;
     a3 = Aql*gama;
     a4 = AquhatI*uhat_analy;
     tetstest=Aqu*u_analy+Aqq*q_analy+Aql*gama-AquhatI*uhat_analy;
%     
     disp(['Cut Element: TesttestElem' num2str(iElem) ]);
     disp(max(abs(tetstest)));
   
    
    %a=AquhatI*uhat_analy
    
% %     
% %     %smalltest
% %     
     u_analy=analiticalSolutionLaplace(Xe,time);
     x = Xe(:,1);
     y = Xe(:,2);
     
     dxu=calculategradx(Xe);
     dyu=calculategrady(Xe);
% %     
     q_analy=zeros(2*length(u_analy),1);
     du_dx = 0+y.*y;
     du_dy = 0+2.*x.*y;
     q_analy(1:2:end)=-1.*du_dx;
     q_analy(2:2:end)=-1.*du_dy;
% %     
% %     
     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama_x=Xe(ce,1);
     coord_gama_y=Xe(ce,2);
     coord_gama=[coord_gama_x,coord_gama_y];
     gama=analiticalSolutionLaplace(coord_gama,time);
% %     
% %     
     test_result_1=fe_s+fe_uns-(Auu*u_analy)-(Auq*q_analy)-(Aul*gama)+(AuuhatI*uhat_analy) ;
     test_result_2= (Aqq*q_analy)+(Aqu*u_analy)+(Aql*gama)-(AquhatI*uhat_analy);
% %     
     disp(['Cut element: Test result equation 1 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_1)));
     disp(['Cut element: Test result equation 2 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_2)));
% %     
% % 
%} 
    
elseif ~isempty(c)
    
    % Elemental mapping
    Auu = Auu_vol_t+Auu_vol_s+Auu_surf;
    Aqu = -Auq'; %Aul = -Alu'; 
    Aql = Alq';
    A = [Auu Auq; Aqu Aqq];
    UQ = -A\[Aul;Aql];
    fUQ= A\[(fe_s+fe_uns);zeros(2*nOfElementNodes,1)];
    U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
    Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q

    %{
% %     %small test
     u_analy=analiticalSolutionLaplace(Xe,time);
     x = Xe(:,1);
     y = Xe(:,2);
% %     
     %dxu=calculategradx(Xe);
     %dyu=calculategrady(Xe);
% %     
     q_analy=zeros(2*length(u_analy),1);
     du_dx = 0+y.*y;
     du_dy = 0+2.*x.*y;
     q_analy(1:2:end)=-1.*du_dx;
     q_analy(2:2:end)=-1.*du_dy;
% %     
     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama_x=Xe(ce,1);
     coord_gama_y=Xe(ce,2);
     coord_gama=[coord_gama_x,coord_gama_y];
     gama=analiticalSolutionLaplace(coord_gama,time);
% %     
% %     
% %     
     test_result_1=fe_s+fe_uns-(Auu*u_analy)-(Auq*q_analy)-(Aul*gama);
     test_result_2=(Aqq*q_analy)+(Aqu*u_analy)+(Aql*gama);
% %     
     disp(['Standard (UnCut) Element: Test result equation 1 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_1)));
     disp(['Standard (UnCut) Element: Test result equation 2 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_2)));
     %}
  
elseif ~isempty(d)
    
    U =zeros(nOfElementNodes,3*nOfFaceNodes);
    Q =zeros(2*nOfElementNodes,3*nOfFaceNodes);
    Uf=zeros(nOfElementNodes,1);
    Qf=zeros(2*nOfElementNodes,1);
    
end












