function [KK,f, QQ, UU, Qf, Uf] = hdg_matrix_PressCorrNeumann(X,T,F,referenceElement,infoFaces,tauP,LS,Elements,dt,Re,a_parm,u,uhat,time)
% This routine is valid only for triangles

nOfElements = size(T,1);
nOfElementNodes = size(referenceElement.NodesCoord,1);
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
    LSe = LS(Te);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)
    %mu=mu_vector(iElem);
    a_cnt = (iElem-1)*nOfElementNodes+(1:nOfElementNodes);
    %
    u_Elem = u(a_cnt,:);
    %
    aux = (1:2*nOfFaceNodes);
    ind = [(Fe(1)-1)*2*nOfFaceNodes + aux,(Fe(2)-1)*2*nOfFaceNodes + aux,(Fe(3)-1)*2*nOfFaceNodes + aux];
    uhat_Elem = uhat(ind);
    
    % elemental matrices
    [Qe,Ue,Qfe,Ufe,Alq,Alp,All] = KKeElementalMatrices(LSe,Elements,Fe,Xe,referenceElement,nOfInteriorFaces,infoFaces,u_Elem,uhat_Elem,tauP,dt,a_parm,iElem,time,Re);
    
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
    Alp=Alp(indL,:);
    All=All(indL,indL);
    
    %The local problem solver is stored for postprocess
    QQ{iElem} = sparse(Qe);  UU{iElem} = sparse(Ue);
    Qf{iElem} = sparse(Qfe); Uf{iElem} = sparse(Ufe);
    
    %Elemental matrices to be assembled
    KKe = Alq*Qe + Alp*Ue + All;
    ffe = -(Alq*Qfe + Alp*Ufe);
    
    %All faces are assembled (Dirichlet included)
    aux = (1:nOfFaceNodes);
    ind = [(Fe(1)-1)*nOfFaceNodes + aux,(Fe(2)-1)*nOfFaceNodes + aux,(Fe(3)-1)*nOfFaceNodes + aux];
    KK(ind,ind)=KK(ind,ind)+ KKe;
    f(ind) = f(ind) + ffe;
end



% ELEMENTAL MATRIX
function [Q,U,Qf,Uf,Alq,Alp,All] = KKeElementalMatrices(LSe,Elements,Fe,Xe,referenceElement,nOfInteriorFaces,infoFaces,u_Elem,uhat_Elem,tau,dt,a_parm,iElem,time,Re)


c=find(iElem==Elements.D1);
d=find(iElem==Elements.D2);


nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)
%tau = [tau,zeros(1,nOfFaces-1)]; %Non-null only for 1st face in each element
%{
l_d = 1;

fe = zeros(nOfElementNodes,1);
Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
Auq = zeros(nOfElementNodes,2*nOfElementNodes);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
Auu = zeros(nOfElementNodes,nOfElementNodes);
Auu_vol = zeros(nOfElementNodes,nOfElementNodes);
Auu_surf = zeros(nOfElementNodes,nOfElementNodes);
Alu = zeros(3*nOfFaceNodes,nOfElementNodes);
Aul = zeros(nOfElementNodes,3*nOfFaceNodes);
All = zeros(3*nOfFaceNodes,3*nOfFaceNodes);
NN_g = zeros(2,2*nOfElementNodes); %Shape functions for vector variable q
NN_xy = zeros(1,2*nOfElementNodes);
%}
%
uhat_Elem_X = zeros(nOfElementNodes,1);
uhat_Elem_Y = zeros(nOfElementNodes,1);
fe = zeros(nOfElementNodes,1);
fe_vol = zeros(nOfElementNodes,1);
fe_surf_1 = zeros(nOfElementNodes,1);
fe_surf_2 = zeros(nOfElementNodes,1);
fe_I = zeros(nOfElementNodes,1);
Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
Apq = zeros(nOfElementNodes,2*nOfElementNodes);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
App = zeros(nOfElementNodes,nOfElementNodes);
%App_surf = zeros(nOfElementNodes,nOfElementNodes);
Alp = zeros(3*nOfFaceNodes,nOfElementNodes);
Apl = zeros(nOfElementNodes,3*nOfFaceNodes);
All = zeros(3*nOfFaceNodes,3*nOfFaceNodes);
NN_g = zeros(2,2*nOfElementNodes); %Shape functions for vector variable q
NN_xy = zeros(1,2*nOfElementNodes);

% Information of the reference element
N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;

if ~isempty(c)   % Standard HDG Element
    
    N = referenceElement.N;
    Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
    IPw = referenceElement.IPweights;
    ngauss=size(IPw,1);
    weight=IPw;
    PtsInt=[];
    FaceInfo=[];
    
    
elseif ~isempty(d)  % Void Element
    
    N=0;
    ngauss=0;
    PtsInt=[];
    FaceInfo=[];
    
else isempty(c) && isempty(d);  % Cut Element
    
    
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

%%% analytical laplacian
%%% sourceTerm = analiticalLaplacianLaplace(Xg,time);

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
    
    %Contribution of the current integration point to the elemental matrix
    Aqq = Aqq + NN_g'*NN_g*dvolu;
    Apq = Apq - N_g'*NN_xy*dvolu; % for N-S eqns.
    %%% Apq = Apq + N_g'*NN_xy*dvolu; % for Laplace eqns.
    %Aqp=  Aqp + NN_g(1,:)'*N_g*dvolu;
    %
    grad_testfunc = [Nx_g ; Ny_g];
    grad_testfunc_T = grad_testfunc';
    %
    u_gauss = N_g*u_Elem;
    %
    fe_vol = fe_vol +grad_testfunc_T*u_gauss'*dvolu/(a_parm*dt); %compute first source term
    %
    %*****OLD*****
    %
    % computing of Auu_vol.
    %c_vec = [c_x c_y];
    %cu = c_vec'*N_g; %cu_t = cu';
    %grad_w = [Nx_g ; Ny_g];
    %grad_w_t = grad_w';
    %Auu_vol_temp = -grad_w_t*cu;
    %App_vol = App_vol + Auu_vol_temp*dvolu; %

    %fe = fe - N_g'*sourceTerm(g)*dvolu;
    %fe = fe + N_g'*sourceTerm(g)*dvolu;
    %
    %sourceTerm = analiticalLaplacianLaplace(Xg(g,:));
    %fe = fe + N_g'*sourceTerm*dvolu; % for Laplace eqns.
    %******************
end


%% Faces computations:
%
% separate x and y components of u_hat velocity.
%
for iface = 1:nOfFaces
    %
    % local numbering of node points on mesh skeleton faces.
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    ind_face1 = (iface-1)*2*nOfFaceNodes + (1:nOfFaceNodes); %% for x comp.
    ind_face2 = ind_face1 + nOfFaceNodes;  %% for y comp.
    %
    uhat_Elem_X(ind_face) = uhat_Elem(ind_face1);
    uhat_Elem_Y(ind_face) = uhat_Elem(ind_face2);
    %
end
%
%% we need to perform the face-flip in order to maintain
%% a coherence between u and u_hat....
%% so that u=u_hat (approximately) for given face...
%% and the essential B.C on the element is satisfied.
%
% Interior faces seen from the second element are flipped to have
% proper orientation
isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)
flipFace = zeros(1,3); %Boolean (1=face to be flipped)
flipFace(isFeInterior)=any(infoFaces.intFaces(Fe(isFeInterior),[1 3])<iElem,2);
indL=1:3*nOfFaceNodes; 
aux_0=nOfFaceNodes:-1:1; 
indflip=[aux_0,nOfFaceNodes+aux_0,2*nOfFaceNodes+aux_0];
aux_1=ones(1,nOfFaceNodes); 
aux_2 = [flipFace(1)*aux_1, flipFace(2)*aux_1, flipFace(3)*aux_1];
indL(aux_2==1)=indflip(aux_2==1); %permutation for local numbering
%
%
uhat_Elem_XN = uhat_Elem_X(indL);
uhat_Elem_YN = uhat_Elem_Y(indL);
%
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces
    
    
    if ~isempty(FaceInfo)
        k=find(FaceInfo(:,1)==iface);
    else
        k=[];
    end
    
    tau_f = tau(iface);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    xf = xe(nodes);
    yf = ye(nodes);
    
    % Inizialization
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    Alq_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    Alq_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    %Alp_f = zeros(nOfFaceNodes,nOfFaceNodes);
    Apl_f = zeros(nOfFaceNodes,nOfFaceNodes);
    App_f = zeros(nOfFaceNodes,nOfFaceNodes);
    All_f = zeros(nOfFaceNodes,nOfFaceNodes);
    fe_uhat_fx = zeros(nOfFaceNodes,1);
    fe_uhat_fy = zeros(nOfFaceNodes,1);
    fe_u_fx = zeros(nOfFaceNodes,1);
    fe_u_fy = zeros(nOfFaceNodes,1);   
    
    
    
    if isempty(k) && any(LSe(faceNodes(iface,:))>0) ;  % NORMAL FACES
        
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
        
        % Shape functions and derivatives at the current Gauss point
        Nf_g2 = Nf_g(g,:);
        Nfxi_g2 = Nfxi_g(g,:);
        
        % Integration weight
        xyDer_g = Nfxi_g2*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline=IPwFace(g)*xyDerNorm_g;
        
        % Unit normal to the boundary
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) -t_g(1)];
        
        %{
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
            App_f = App_f + tau_f*Nf_g(g,:)'*Nf_g(g,:)*dline; % Auu_surf
            Apl_f = Apl_f + bb*Nf_g(g,:)'*Nf_g(g,:)*dline; 
            All_f = All_f + bb*Nf_g(g,:)'*Nf_g(g,:)*dline;
        else
            disp('Lemma 3.1 violated');
        end  
        %}

        % Contribution of the current integration point to the elemental matrix
        %Alq_fnx = Alq_fnx + Nf_g2'*Nf_g2*n_g(1)*dline;
        %Alq_fny = Alq_fny + Nf_g2'*Nf_g2*n_g(2)*dline;
        %App_f = App_f + tau_f *Nf_g2'*Nf_g2*dline;
        %All_f = All_f - tau_f *Nf_g2'*Nf_g2*dline;
        %
        %
        % Contribution of the current Gauss point to the elemental matrices
        %   
        %%%All_f = All_f - tau_f*(Nf_g2')*Nf_g2*dline; % for Laplace eqns.
        All_f = All_f + tau_f*(Nf_g2')*Nf_g2*dline; % for N-S eqns.
        %
        Alq_fnx = Alq_fnx + Nf_g2'*Nf_g2*n_g(1)*dline;
        Alq_fny = Alq_fny + Nf_g2'*Nf_g2*n_g(2)*dline;
        %
        %Alu_f = Alu_f + tau_f*(Nf_g')*Nf_g*dline; 
        %
        App_f = App_f + tau_f*(Nf_g2')*Nf_g2*dline;
        %
        Apl_f = Apl_f - tau_f*(Nf_g2)'*Nf_g2*dline;
        %
        %compute the second source term
        %uhatX_gaussFace = Nf_g*uhat_Elem(ind_face1);
        %uhatY_gaussFace = Nf_g*uhat_Elem(ind_face2); 
        %
        uhatX_gaussFace = Nf_g2*uhat_Elem_XN(ind_face);
        uhatY_gaussFace = Nf_g2*uhat_Elem_YN(ind_face);        
        %
        dd_hat_x = (Nf_g2'*n_g(1))*uhatX_gaussFace*dline/(a_parm*dt);
        dd_hat_y = (Nf_g2'*n_g(2))*uhatY_gaussFace*dline/(a_parm*dt);
        %
        fe_uhat_fx = fe_uhat_fx - dd_hat_x;
        fe_uhat_fy = fe_uhat_fy - dd_hat_y;
        %
        %
        uX_gaussFace = Nf_g2*u_Elem(nodes,1);
        uY_gaussFace = Nf_g2*u_Elem(nodes,2);
        %
        dd_x = (Nf_g2'*n_g(1))*uX_gaussFace*dline/(a_parm*dt);
        dd_y = (Nf_g2'*n_g(2))*uY_gaussFace*dline/(a_parm*dt);
        %
        % 
        fe_u_fx = fe_u_fx + dd_x;
        fe_u_fy = fe_u_fy + dd_y;
        %                  
    end
    
    %{
    Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Alq_fnx;
    Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Alq_fny;
    App_surf(nodes,nodes) = App_surf(nodes,nodes) + App_f; %
    Apl(nodes,ind_face) = Apl(nodes,ind_face) + Apl_f; %
    Alp(ind_face,nodes) = Alp(ind_face,nodes) + App_f;
    All(ind_face,ind_face) = All(ind_face,ind_face) + All_f;
    %}
    %
    Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Alq_fnx;
    Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Alq_fny;
    %
    %Alu(ind_face,nodes) = Alu(ind_face,nodes) + Alu_f;
    %
    Apl(nodes,ind_face) = Apl(nodes,ind_face) + Apl_f;
    %
    App(nodes,nodes) = App(nodes,nodes) + App_f;
    %
    Alp(ind_face,nodes) = Alp(ind_face,nodes) + App_f;
    %
    All(ind_face,ind_face) = All(ind_face,ind_face) + All_f;
    %
    fe_surf_1(nodes) = fe_surf_1(nodes) + (fe_uhat_fx+fe_uhat_fy) + (fe_u_fx+fe_u_fy); 
    %    
    fe_surf_2(nodes) = fe_surf_2(nodes) - (fe_u_fx+fe_u_fy); 
    %    
end


if  isempty(c) && isempty (d) % for interface I of the cut element
    
    %Calculation of new AquuhatI,AuuhatI and AuuI matrices because of applied
    %Neumann boundary conditions on I
    
    % AAAANNNNDDDD
       
    %Calculation of new AuhatqI,AuhatuI and AuhatuhatI matrices because of applied
    %Neumann boundary conditions on I
    
    
    p=referenceElement.degree;
    g=length(referenceElement.IPweights1d);
    AqphatI=zeros(2*nOfElementNodes,nOfFaceNodes);
    
    
    PtsInt=reshape(PtsInt',2,p+1)'; %interface p+1 nodes in the REFERENCE element
    
    zg=referenceElement.N1d*PtsInt;%integration points on the interface in the REFERENCE element

    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zg);
    Nzg = shapeFunctions(:,:,1)';  %2d shape functions for the Gauss points on interface I %%(on the reference element at integration points)
    
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,PtsInt);
    NPtsInt = shapeFunctions(:,:,1)';  %2D Shape functions on the REFERENCE element at interface nodes
    
    Pphy=NPtsInt*Xe; %p+1 interface nodes on the REAL element
    Pgauss = referenceElement.N1d*Pphy;

    
%     figure(1); hold on,
%     
%     plot(Pphy(:,1),Pphy(:,2),'ro')
    
    Iprime=referenceElement.N1dxi*Pphy;
    I=referenceElement.N1d*Pphy;
    
%     figure(1); hold on,
%     
%     quiver(I(:,1),I(:,2),Iprime(:,1),Iprime(:,2))
    
    %Initialization
    
    ApphatI=zeros(nOfElementNodes,nOfFaceNodes);
    %AphatpI=zeros(nOfFaceNodes,nOfElementNodes);
    AppI=zeros(nOfElementNodes,nOfElementNodes);
    %AqphatI2 = zeros(2*nOfElementNodes,nOfFaceNodes); 
    AqphatI_x_f=zeros(nOfElementNodes,nOfFaceNodes);
    AqphatI_y_f=zeros(nOfElementNodes,nOfFaceNodes);
    AphatphatI=zeros(nOfFaceNodes,nOfFaceNodes);
    ge=zeros(nOfFaceNodes,1);
    %ux_hatI = zeros(g,1);
    %uy_hatI = zeros(g,1);    
    %
    for igauss=1:g
        
        
        % Integration weight
        normIprime = norm(Iprime(igauss,:));
        wg=referenceElement.IPweights1d(igauss)*normIprime;
        
        % Unit normal to the boundary
        t_g = Iprime(igauss,:)/normIprime;
        n_g = [t_g(2) -t_g(1)];
        %
        %
        %ApphatI=ApphatI+bb*Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*wg;
        gaga = referenceElement.N1d(igauss,:);
        haha = Nzg(igauss,:);
        fafa = 0.*haha;
        %
        cg_5I_ng_temp = [haha fafa; fafa haha];
        cg_5I_ng = n_g*cg_5I_ng_temp;
        %d = b*c;
        %ApphatI=ApphatI-tau_f*Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*wg;
        ApphatI=ApphatI-tau_f*haha'*gaga*wg; 
        %AphatpI=AphatpI-tau_f*referenceElement.N1d(igauss,:)'*Nzg(igauss,:)*wg;
        %
        %AppI=AppI+tau_f*Nzg(igauss,:)'*Nzg(igauss,:)*wg;
        AppI=AppI+tau_f*haha'*haha*wg;
        %
        %% for N-S eqns.
        %AqphatI_x_f=AqphatI_x_f-Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(1)*wg;
        %AqphatI_y_f=AqphatI_y_f-Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(2)*wg;
        AqphatI_x_f=AqphatI_x_f-haha'*gaga*n_g(1)*wg;
        AqphatI_y_f=AqphatI_y_f-haha'*gaga*n_g(2)*wg;
        %
        %% for Laplace eqns.
        %%%AqphatI_x_f=AqphatI_x_f+Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(1)*wg;
        %%%AqphatI_y_f=AqphatI_y_f+Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(2)*wg;
        %% for N-S eqns.
        %AphatphatI=AphatphatI+tau_f*referenceElement.N1d(igauss,:)'*referenceElement.N1d(igauss,:)*wg;
        AphatphatI = AphatphatI + tau_f*gaga'*gaga*wg;
        %% for Laplace eqns.
        %%%AphatphatI=AphatphatI-tau_f*referenceElement.N1d(igauss,:)'*referenceElement.N1d(igauss,:)*wg;
        %ge= ge+referenceElement.N1d(igauss,:)'*-NeumannCondition(n_g)*wg;
        %
        %aa = referenceElement.N1d(igauss,:)';
        bb = NeumannConditionPcorr(Pgauss(igauss,:),n_g);
        ge = ge + gaga'*bb*wg;
        %
        %%%ge=ge+referenceElement.N1d(igauss,:)'*-NeumannCondition(n_g,Pgauss(igauss,:),c_x,c_y)*wg;
        %        
              
    %end
    %
        %
    %for igauss=1:g
        %       
        % Integration weight
        %normIprime = norm(Iprime(igauss,:));
        %wg=referenceElement.IPweights1d(igauss)*normIprime;
        %
        % Unit normal to the boundary
        %t_g = Iprime(igauss,:)/normIprime;
        %n_g = [t_g(2) -t_g(1)];
        %% CHECK IT AS WELL
        %
        uhatI_g = analyticalVelocityStokes(Pgauss(igauss,:),time,Re);
        %gI_dot_ng = n_g*uhatI_g;
        %%%fe_I_temp = fe_I_temp -referenceElement.N1d(igauss,:)'*(n_g*uhatI_g)*wg/(a_parm*dt);
        %fe_I = fe_I -Nzg(igauss,:)'*(n_g*uhatI_g)*wg/(a_parm*dt);
        fe_I = fe_I - haha'*(n_g*uhatI_g)*wg/(a_parm*dt);
        %      
        %     
    end
    %
    %{
    %*****
    dline = spdiags(dxdxiNorm.*IPw_f,0,g,g);
    M = (referenceElement.N1d)'*(dline*referenceElement.N1d);
    b1 = (referenceElement.N1d)'*(dline*ux_hatI);
    b2 = (referenceElement.N1d)'*(dline*uy_hatI);
    uf_x = M\b1;
    uf_y = M\b2;
    uface = [uf_x ; uf_y];
    %*****
    %}
    %
    AqphatI(1:2:end,:) = AqphatI_x_f;
    AqphatI(2:2:end,:) = AqphatI_y_f;
     
    AphatqI=-AqphatI'; % for N-S eqns.
    %%%AphatqI=AqphatI'; % for Laplace eqns.
    %
    AphatpI=ApphatI'; % for N-S eqns.
    %%%AphatpI=-ApphatI'; % for Laplace eqns.
    
    Me=(AphatphatI)^-1;
    M2_e = inv(AphatphatI);
    Tq=-Me*AphatqI;
    Tp=-Me*AphatpI;
    
    %te=-Me*ge;
    te=Me*ge;
    %
    Aqp =-Apq';
    Aql = -Alq'; 
    Alp = Apl'; 
   
    Appnew=App+AppI+ApphatI*Tp;
    Apqnew=Apq+ApphatI*Tq;
    %
    Aqpnew=Aqp+AqphatI*Tp;
    Aqqnew=Aqq+AqphatI*Tq;
    
    %Elemental mapping
    %
    fe = fe_vol + fe_surf_1 + fe_surf_2 + fe_I - (ApphatI*te);
    fq = -1.*(AqphatI*te);
    %fq = zeros(2*nOfElementNodes,1);
    %
    %Aqp = -Apq';    
    %Aqu = -Auq';
    %Aul = -Alu'; 
    %Aql = Alq';
    
    A = [Appnew Apqnew; Aqpnew Aqqnew];
    UQ = -A\[Apl;Aql];
    %%% fUQ= A\[(fe-(ApphatI*te));-(AqphatI*te)]; % for Laplace eqn.
    fUQ= A\[fe;fq];
    U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
    Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q
    
     %smalltest1
     %{
     u_analy=analiticalSolutionLaplace(Xe,time);
     x = Xe(:,1);
     y = Xe(:,2);
     du_dx = 0+1.*(y.^3);
     du_dy = 0+3.*x.*(y.^2);
     %***compute du/dx***
     %trmA = sin(pi*x)+pi.*cos(pi*x);
     %du_dx = exp(x+y).*sin(pi*y).*trmA;
     %***compute du/dy***
     % trmB = sin(pi*y)+pi.*cos(pi*y);
     % du_dy = exp(x+y).*sin(pi*x).*trmB;
      %*******************
     q_analy=zeros(2*length(u_analy),1);
     q_analy(1:2:end)=1.*du_dx;
     q_analy(2:2:end)=1.*du_dy;
    
% % %     
     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama_x=Xe(ce,1);
     coord_gama_y=Xe(ce,2);
     coord_gama=[coord_gama_x,coord_gama_y];
     gama=analiticalSolutionLaplace(coord_gama,time);
% % %     
% % %     
% % % 
     a1 = Appnew*u_analy;
     a2 = Apqnew*q_analy;
     a3 = Apl*gama;
     a4 = ApphatI*te;    
     test_result_1=fe-(Appnew*u_analy)-(Apqnew*q_analy)-(Apl*gama)-(ApphatI*te) ;
     %test_result_1=fe-(Appnew*u_analy)-(Apqnew*q_analy)-(Apl*gama);

     b1 = Aqqnew*q_analy;
     b2 = Aqpnew*u_analy;
     b3 = Aql*gama;
     b4 = AqphatI*te;
     test_result_2= (Aqqnew*q_analy)+(Aqpnew*u_analy)+(Aql*gama)+(AqphatI*te);
     %test_result_2= (Aqqnew*q_analy)+(Aqpnew*u_analy)+(Aql*gama);

% % %     
     disp(['Test result equation 1 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_1)));
     disp(['Test result equation 2 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_2)));
% % %     
% % %     %smalltest2
% % %      
     uhat_analy= analiticalSolutionLaplace(Pphy,time);
     c1 = AphatqI*q_analy;
     c2 = AphatpI*u_analy;
     c3 = AphatphatI*uhat_analy;
     c4 = ge;
     test=c1+c2+c3-ge;
% % %     
     disp(['NewTest' num2str(iElem) ]);
     disp(test);
% % %     
% % %     %smalltest3
% % %     
% % %     
     d1 = Tq*q_analy;
     d2 = Tp*u_analy;
     %d3 = AuhatuhatI*uhat_analy;

     %test3=d1+d2-te-uhat_analy;
     test3=d1+d2+te-uhat_analy;
     %test3=d1+d2-uhat_analy;
     
     disp(['NewNewTest' num2str(iElem) ]);
     disp(test3);
     disp("Hola")
     %}

elseif ~isempty(c) % for standard HDG Element
    
    % Elemental mapping
    %
    Aqp = -Apq';
    Aql = -Alq'; 
    Alp = Apl'; 
    %   
    %
    fe = fe_vol + fe_surf_1 + fe_surf_2;
    fq = zeros(2*nOfElementNodes,1);
    %
    A = [App Apq; Aqp Aqq];
    UQ = -A\[Apl;Aql];
    fUQ= A\[fe;fq];
    U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
    Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q

    
%     %small test
%{
     u_analy=analiticalSolutionLaplace(Xe,time);
     x = Xe(:,1);
     y = Xe(:,2);
     du_dx = 0+1.*y.^3;
     du_dy = 0+3*x.*y.^2;
     q_analy=zeros(2*length(u_analy),1);
     q_analy(1:2:end)=1.*du_dx;
     q_analy(2:2:end)=1.*du_dy;
%     
     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama_x=Xe(ce,1);
     coord_gama_y=Xe(ce,2);
     coord_gama=[coord_gama_x,coord_gama_y];
     gama=analiticalSolutionLaplace(coord_gama,time);
%     
%     
%     
     test_result_1=fe-(App*u_analy)-(Apq*q_analy)-(Apl*gama);
     test_result_2=(Aqq*q_analy)+(Aqp*u_analy)+(Aql*gama);
%     
     disp(['Test result equation 1 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_1)));
     disp(['Test result equation 2 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_2)));
     disp("Hola")
%}
    
    
elseif ~isempty(d) % element in the Void region
    
    U =zeros(nOfElementNodes,3*nOfFaceNodes);
    Q =zeros(2*nOfElementNodes,3*nOfFaceNodes);
    Uf=zeros(nOfElementNodes,1);
    Qf=zeros(2*nOfElementNodes,1);
    
end

%disp('Hola')












