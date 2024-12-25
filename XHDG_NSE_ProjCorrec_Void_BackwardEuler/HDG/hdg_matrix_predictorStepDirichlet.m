
function [KK,f, QQ, UU, Qf, Uf] = hdg_matrix_predictorStepDirichlet(X,T,F,referenceElement,infoFaces,tau,LS,Elements,dt,cnt,u,p,time,Re,a_parm)
% This routine is valid only for triangles

nOfElements = size(T,1);
nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(infoFaces.intFaces,1)+size(infoFaces.extFaces,1);
%
%nDOF = nOfFaces*nOfFaceNodes;
%KK = spalloc(nDOF,nDOF,nOfFaceNodes*nDOF);
%f = zeros(nDOF,1);
%
nDOF = 2*nOfFaces*nOfFaceNodes;
KK = spalloc(nDOF,nDOF,2*nOfFaceNodes*nDOF);
f = zeros(nDOF,1);
%
QQ = cell(nOfElements,1); UU = cell(nOfElements,1); Qf = cell(nOfElements,1); Uf = cell(nOfElements,1);
%
aux=nOfFaceNodes:-1:1; aux=[aux,aux+nOfFaceNodes]; 
indflip=[aux,2*nOfFaceNodes+aux,4*nOfFaceNodes+aux];

% loop in elements

for iElem = 1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    LSe = LS(Te);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)
    %u_Elem = u((iElem-1)*nOfElementNodes+(1:nOfElementNodes));
    a_cnt = (iElem-1)*nOfElementNodes+(1:nOfElementNodes);
    %
    u_intm = u(a_cnt,:);
    %
    u_Elem = [u_intm(:,1) ; u_intm(:,2)];
    p_Elem = p(a_cnt);
    %u_Elem = [u(a_cnt,1) ; u(a_cnt,2)];
    
    % elemental matrices
    [Qe,Ue,Qfe,Ufe,Alq,Alu,All,gn_press] = KKeElementalMatrices(Xe,referenceElement,LSe,Elements,tau,dt,cnt,u_Elem,p_Elem,time,Re,a_parm,iElem);
    %{
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
    %}
    % Collect u_I for all the Elements.
    %uhat_I(iElem,1) = u_I; 
    %
    % Interior faces seen from the second element are flipped to have
    % proper orientation
    flipFace = zeros(1,3); %Boolean (1=face to be flipped)
    flipFace(isFeInterior)=any(infoFaces.intFaces(Fe(isFeInterior),[1 3])<iElem,2);
    indL=1:(2*3*nOfFaceNodes); 
    %aux=nOfFaceNodes:-1:1; indflip=[aux,nOfFaceNodes+aux,2*nOfFaceNodes+aux];
    aux=ones(1,2*nOfFaceNodes); aux = [flipFace(1)*aux, flipFace(2)*aux, flipFace(3)*aux];
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
    ffe = gn_press-(Alq*Qfe + Alu*Ufe);
    
    %All faces are assembled (Dirichlet included)
    aux = (1:2*nOfFaceNodes);
    ind = [(Fe(1)-1)*2*nOfFaceNodes + aux,(Fe(2)-1)*2*nOfFaceNodes + aux,(Fe(3)-1)*2*nOfFaceNodes + aux];
    KK(ind,ind)=KK(ind,ind)+ KKe;
    f(ind) = f(ind) + ffe;
    %
    %disp(iElem)
    %    
end



% ELEMENTAL MATRIX
function [Q,U,Qf,Uf,Alq,Alu,All,gn_press] = KKeElementalMatrices(Xe,referenceElement,LSe,Elements,tau,dt,cnt,u_Elem,p_Elem,time,Re,a_parm,iElem)


c=find(iElem==Elements.D1);
d=find(iElem==Elements.D2);
%{
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
%}

% Information of the reference element
%N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;

nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)
%% *****OLD*****
%tau = [tau,zeros(1,nOfFaces-1)]; %Non-null only for 1st face in each element
%N = referenceElement.N;
%Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
%IPw = referenceElement.IPweights; ngauss=size(IPw,1);
%Xg = N*Xe; %physical coordinates
%% ************

% Information of the reference element
N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;

% x and y coordinates of the element nodes
%xe = Xe(:,1); ye = Xe(:,2);
%
aux1=1:nOfElementNodes; aux2=aux1+nOfElementNodes; aux3=aux2+nOfElementNodes; aux4=aux3+nOfElementNodes;

%%%***Allocate the arrays Aii used in vol. integrals of local elemental problem***
Auq = zeros(2*nOfElementNodes,4*nOfElementNodes);
Aqq = zeros(4*nOfElementNodes,4*nOfElementNodes);
Auu_vol_t = zeros(2*nOfElementNodes,2*nOfElementNodes);
fe_s_vel = zeros(2*nOfElementNodes,1);
fe_s_press_vol = zeros(2*nOfElementNodes,1);
fe_s_press_surf = zeros(2*nOfElementNodes,1);
fe_uns = zeros(2*nOfElementNodes,1);
%Shape functions for vector variable q
NN_g = zeros(4,4*nOfElementNodes); 
NN_xy = zeros(1,2*nOfElementNodes);
NNN_xy = zeros(2,4*nOfElementNodes);


if ~isempty(c) % Standard HDG element
    
    N = referenceElement.N;
    Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
    IPw = referenceElement.IPweights;
    ngauss=size(IPw,1);
    weight=IPw;
    PtsInt=[];
    FaceInfo=[];
        
elseif ~isempty(d) % Void Element
    
    N=0;
    ngauss=0;
    PtsInt=[];
    FaceInfo=[];
    
else isempty(c) && isempty(d); % Cut Element
    
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
%Xg = N*Xe; %physical coordinates
IPw_f = referenceElement.IPweights1d;

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

%% ********************************************************
%% OLD*****
%% ********************************************************
% analytical laplacian
%sourceTerm = analiticalLaplacianLaplace(Xg,c_x,c_y,time);
%{
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
%}
%% ********************************************************
%% ********************************************************
%
%% Volume computation: loop in Gauss points
for g = 1:ngauss
    
    %Shape functions and derivatives at the current Gauss point
    N_g = N(g,:);
    %NN_g(1,1:2:end)=N_g; NN_g(2,2:2:end)=N_g;
    NN_g(1,aux1)=N_g; NN_g(2,aux2)=N_g; NN_g(3,aux3)=N_g; NN_g(4,aux4)=N_g;
    Nxi_g = Nxi(g,:);    Neta_g = Neta(g,:);
    
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end

    %Integration weight
    dvolu=weight(g)*det(J);
    
    %physical coordinates of the current Gauss point on the Element
    xy_g=N_g*Xe;
 
    %x and y derivatives
    invJ = inv(J);  
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    
    %NN_xy(1:2:end) = Nx_g;    NN_xy(2:2:end) = Ny_g;
    NN_xy(aux1) = Nx_g;    NN_xy(aux2) = Ny_g;
    
    NNN_xy(1,aux1)=Nx_g; NNN_xy(1,aux2)=Ny_g; NNN_xy(2,aux3)=Nx_g; NNN_xy(2,aux4)=Ny_g;
    
    cg_1x = [N_g ; 0.*N_g]';
    cg_1y = [0.*N_g ; N_g]';
    cg_1 = [cg_1x ; cg_1y];
    
    cg_2x = [N_g ; 0.*N_g];
    cg_2y = [0.*N_g ; N_g];
    cg_2 = [cg_2x  cg_2y];
     
    %Contribution of the current Gauss point to the elemental matrix
    Aqq = Aqq + NN_g'*NN_g*dvolu*Re;
    Auq = Auq - cg_1*NNN_xy*dvolu;  
    Auu_vol_t = Auu_vol_t + (cg_1*cg_2*dvolu/(a_parm*dt));
     
    cg_0 = sourceTermStokes(xy_g,Re,time);
    fe_s_vel = fe_s_vel + cg_1*cg_0'*dvolu;
    
    if(cnt==1)
        p_old_gauss = press_initial_condition(xy_g,time);
    else
        p_old_gauss = N_g*p_Elem;
    end 
    %
    %
    fe_s_press_vol = fe_s_press_vol + NN_xy'*p_old_gauss*dvolu; %first pressure source term
    
    
    if(cnt==1)
        u_old_gauss = vel_initial_condition(xy_g);
    else
        u_old_gauss = cg_2*u_Elem; 
    end 
    %
    %fe_uns = fe_uns + N_g'*u_old_gauss*dvolu/dt;
    fe_uns = fe_uns + cg_1*u_old_gauss*dvolu/(a_parm*dt);
    
end
%inv_Aqq = inv(Aqq);
% Elemental mapping for terms containing volume integrals
Aqu = -Auq';

%%%***Allocate the arrays Aii used in surf. integrals of local elemental problem***
Auu_surf = zeros(2*nOfElementNodes,2*nOfElementNodes);
Alq = zeros(2*3*nOfFaceNodes,4*nOfElementNodes);
Aql = zeros(4*nOfElementNodes,2*3*nOfFaceNodes);
Alu = zeros(2*3*nOfFaceNodes,2*nOfElementNodes);
Aul = zeros(2*nOfElementNodes,2*3*nOfFaceNodes);
All = zeros(2*3*nOfFaceNodes,2*3*nOfFaceNodes);
gn_press = zeros(2*3*nOfFaceNodes,1);

%% Faces computations:
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces
    
    
    if ~isempty(FaceInfo)
        k=find(FaceInfo(:,1)==iface);
    else
        k=[];
    end
    
    tau_f = tau(iface);
    
    % Nodes in the face
    %nodes = faceNodes(iface,:);
    %xf = xe(nodes);
    %yf = ye(nodes);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    nodes1 = faceNodes(iface,:);
    nodes2 = nodes1+nOfElementNodes; nodes3=nodes2+nOfElementNodes; nodes4=nodes3+nOfElementNodes;
    
    %ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);  
    %node point numbering for mesh skeleton faces.
    ind_face1 = (iface-1)*2*nOfFaceNodes + (1:nOfFaceNodes); %% for x comp.
    ind_face2 = ind_face1 + nOfFaceNodes;  %% for y comp.
    
    xf = xe(nodes);
    yf = ye(nodes);
       
    % Inizialization
    %ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    %Alq_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    %Alq_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    %Auu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    %Aul_f = zeros(nOfFaceNodes,nOfFaceNodes);
    %All_f = zeros(nOfFaceNodes);
    %allocate local arrays for local computations on the faces.
    Auu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    Aql_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    Aql_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_f = zeros(nOfFaceNodes,nOfFaceNodes);
    All_f = zeros(nOfFaceNodes,nOfFaceNodes);   
    %Alu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    fe_s_press_fx = zeros(nOfFaceNodes,1);
    fe_s_press_fy = zeros(nOfFaceNodes,1);    
    
    
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
        
        % Gauss points position on the current face
        xfg = Nf_g*xf;
        yfg = Nf_g*yf;  
        
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
        
        % physical location of the Gauss points
        xy_gf = [xfg(g) yfg(g)];
        %
        % Contribution of the current integration point to the elemental matrix
        %
         Auu_f = Auu_f + tau_f*(Nf_g2')*Nf_g2*dline;
         Aql_fnx = Aql_fnx - Nf_g2'*Nf_g2*n_g(1)*dline;
         Aql_fny = Aql_fny - Nf_g2'*Nf_g2*n_g(2)*dline;
         Aul_f = Aul_f - tau_f*(Nf_g2')*Nf_g2*dline;
         All_f = All_f + tau_f*(Nf_g2')*Nf_g2*dline;
         
         %compute the second pressure source term 
         if(cnt==1)
             p_old_gaussFace = press_initial_condition(xy_gf,time);
         else
             p_old_gaussFace = Nf_g2*p_Elem(nodes);
         end
         dd_x = (Nf_g2'*n_g(1))*p_old_gaussFace*dline;
         dd_y = (Nf_g2'*n_g(2))*p_old_gaussFace*dline;
         % for the use in local eqn. of predictor velocity (in source term)
         % and for the use in the global eqn. (in source term)
         fe_s_press_fx = fe_s_press_fx + dd_x;
         fe_s_press_fy = fe_s_press_fy + dd_y; 
         %        
        %% *********************************************
        %% OLD *******
        %% *********************************************
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
            Auu_f = Auu_f + tau_f*Nf_g(g,:)'*Nf_g(g,:)*dline; % Auu_surf
            Aul_f = Aul_f + bb*Nf_g(g,:)'*Nf_g(g,:)*dline; 
            All_f = All_f + bb*Nf_g(g,:)'*Nf_g(g,:)*dline;
        else
            disp('Lemma 3.1 violated');
        end  
        %}
        %**************************************************
        %**************************************************
        %
    end
    %{    
    Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Alq_fnx;
    Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Alq_fny;
    Auu_surf(nodes,nodes) = Auu_surf(nodes,nodes) + Auu_f; %
    Aul(nodes,ind_face) = Aul(nodes,ind_face) + Aul_f; %
    Alu(ind_face,nodes) = Alu(ind_face,nodes) + Auu_f;
    All(ind_face,ind_face) = All(ind_face,ind_face) + All_f;
    %}
    %    
    Auu_surf(nodes1,nodes1) = Auu_surf(nodes1,nodes1)+Auu_f; 
    Auu_surf(nodes2,nodes2) = Auu_surf(nodes2,nodes2)+Auu_f; 
    %
    Aul(nodes1,ind_face1) = Aul(nodes1,ind_face1)+Aul_f; 
    Aul(nodes2,ind_face2) = Aul(nodes2,ind_face2)+Aul_f;
    %
    Aql(nodes1,ind_face1) = Aql(nodes1,ind_face1)+Aql_fnx;
    Aql(nodes2,ind_face1) = Aql(nodes2,ind_face1)+Aql_fny;
    Aql(nodes3,ind_face2) = Aql(nodes3,ind_face2)+Aql_fnx;
    Aql(nodes4,ind_face2) = Aql(nodes4,ind_face2)+Aql_fny; 
    %
    All(ind_face1,ind_face1)=All(ind_face1,ind_face1)+All_f;
    All(ind_face2,ind_face2)=All(ind_face2,ind_face2)+All_f;
    %
    fe_s_press_surf(nodes1) = fe_s_press_surf(nodes1)+fe_s_press_fx;
    fe_s_press_surf(nodes2) = fe_s_press_surf(nodes2)+fe_s_press_fy;
    %
    gn_press(ind_face1) = gn_press(ind_face1)+fe_s_press_fx;
    gn_press(ind_face2) = gn_press(ind_face2)+fe_s_press_fy;
    %    
end


if  isempty(c) && isempty (d) % for interface I of the cut element
    
    %Calculation of new AquuhatI,AuuhatI and AuuI matrices because of applied
    %Dirichlet boundary conditions on I
    
    p=referenceElement.degree;
    % no. of Gauss points on the Interface I;
    g=length(referenceElement.IPweights1d);
    %AquhatI=zeros(2*nOfElementNodes,nOfFaceNodes);
    AquhatI=zeros(4*nOfElementNodes,2*nOfFaceNodes);
    %Aql = zeros(4*nOfElementNodes,2*3*nOfFaceNodes);
    
    
    PtsInt=reshape(PtsInt',2,p+1)'; %interface p+1 nodes in the REFERENCE element
    
    N_1d_tmp = referenceElement.N1d; N_1dxi_tmp = referenceElement.N1dxi;
    
    %zg=referenceElement.N1d*PtsInt;%integration points on the interface in the REFERENCE element
    zg=N_1d_tmp*PtsInt;%integration points on the interface in the REFERENCE element
    
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zg);
    Nzg = shapeFunctions(:,:,1)'; %2d shape functions for the Gauss points on interface I %%(on the reference element at integration points)
    
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,PtsInt);
    NPtsInt = shapeFunctions(:,:,1)';  %2D Shape functions on the REFERENCE element at interface nodes
    NPtsInt_1 = NPtsInt';
    
    Pphy=NPtsInt*Xe; % physical nodes on the interface in the cut element
    %Pgauss = referenceElement.N1d*Pphy;
    Pgauss = N_1d_tmp*Pphy; % location of Gauss points on the interface
       
    %Iprime=referenceElement.N1dxi*Pphy;
    Iprime=N_1dxi_tmp*Pphy;
    % I=referenceElement.N1d*Pphy;
    %dxdxi = Nx1d*xf; dydxi = Nx1d*yf;
    dxdxiNorm = sqrt(Iprime(:,1).^2+Iprime(:,2).^2);
    
    %Initialization
    
    %AuuhatI=zeros(nOfElementNodes,nOfFaceNodes);
    %AuuI=zeros(nOfElementNodes,nOfElementNodes);
    %AquhatI_x_f=zeros(nOfElementNodes,nOfFaceNodes);
    %AquhatI_y_f=zeros(nOfElementNodes,nOfFaceNodes);
    %
    %% These matrices are on the interface but they add-up into the
    %% matrices for the whole element, hence their allocation also contains
    %% the nOfElementNodes....
    %
    AuphatI = zeros(2*nOfElementNodes,nOfFaceNodes); 
    AuuhatI=zeros(2*nOfElementNodes,2*nOfFaceNodes);
    AuuI=zeros(2*nOfElementNodes,2*nOfElementNodes);
    AquhatI_x_f=zeros(2*nOfElementNodes,2*nOfFaceNodes);
    AquhatI_y_f=zeros(2*nOfElementNodes,2*nOfFaceNodes);
    ux_hatI = zeros(g,1);
    uy_hatI = zeros(g,1);
    %
    for igauss=1:g
        
        
        % Integration weight
        normIprime = norm(Iprime(igauss,:));
        wg=referenceElement.IPweights1d(igauss)*normIprime;
        
        % Unit normal to the boundary
        t_g = Iprime(igauss,:)/normIprime;
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
        %}

        %{
        AuuhatI = AuuhatI - tau_f*Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*wg;
        AuuI = AuuI + tau_f*Nzg(igauss,:)'*Nzg(igauss,:)*wg;
        AquhatI_x_f =  AquhatI_x_f + Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(1)*wg;
        AquhatI_y_f =  AquhatI_y_f + Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(2)*wg;
        %}
        
        cg_1_Ix = [Nzg(igauss,:) ; 0.*Nzg(igauss,:)]';
        cg_1_Iy = [0.*Nzg(igauss,:) ; Nzg(igauss,:)]';
        cg_1_I = [cg_1_Ix ; cg_1_Iy];
    
        cg_2_Ix = [Nzg(igauss,:) ; 0.*Nzg(igauss,:)];
        cg_2_Iy = [0.*Nzg(igauss,:) ; Nzg(igauss,:)];
        cg_2_I = [cg_2_Ix  cg_2_Iy];  
        
        gaga = N_1d_tmp(igauss,:);
        cg_3_I = [gaga 0.*gaga; 0.*gaga gaga];

        haha = Nzg(igauss,:);
        fafa = 0.*haha;
        cg_4_I1 = [haha fafa fafa fafa; fafa haha fafa fafa];
        cg_4_I2 = [fafa fafa haha fafa; fafa fafa fafa haha];
        %
        cg_4I_ng1 = n_g*cg_4_I1; 
        cg_4I_ng2 = n_g*cg_4_I2;
        %
        cg_4I_ng = [cg_4I_ng1;cg_4I_ng2];

        cg_5I_ng_temp = [haha fafa; fafa haha];
        cg_5I_ng = n_g*cg_5I_ng_temp;
        
        AuuI = AuuI + tau_f*cg_1_I*cg_2_I*wg;
        %%
        AuuhatI = AuuhatI + tau_f*cg_1_I*cg_3_I*wg;
        %%
        AuphatI = AuphatI + cg_5I_ng'*gaga*wg;
        %%
        AquhatI = AquhatI + cg_4I_ng'*cg_3_I*wg;
        %%
        %AquhatI_x_f =  AquhatI_x_f + cg_1_I*cg_3_I*n_g(1)*wg;
        %AquhatI_y_f =  AquhatI_y_f + cg_1_I*cg_3_I*n_g(2)*wg; 
        %
        uhatI_g = analyticalVelocityStokes(Pgauss(igauss,:),time);
        ux_hatI(igauss) = uhatI_g(1); uy_hatI(igauss) = uhatI_g(2);
        %
        %aaa = tau_f*Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*wg;
        %AuuhatI = AuuhatI - tau_f*Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*wg;
        %AuuI = AuuI + tau_f*Nzg(igauss,:)'*Nzg(igauss,:)*wg;
        %AquhatI_x_f =  AquhatI_x_f - Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(1)*wg;
        %AquhatI_y_f =  AquhatI_y_f - Nzg(igauss,:)'*referenceElement.N1d(igauss,:)*n_g(2)*wg;

    end
    
    %AquhatI(1:2:end,:) = AquhatI_x_f;
    %AquhatI(2:2:end,:) = AquhatI_y_f;
    
    %% for Laplace eqns.
    %uhatI_e_old=analiticalSolutionLaplace(Pphy,time);
    %uhatI_e=analiticalSolutionLaplace(Pgauss,time);
    %% for N-S eqns.
    %uhatI_e_old=analyticalVelocityStokes(Pphy,time);
     
    %ux_hatI = uhatI_e(1:2:end);
    %uy_hatI = uhatI_e(2:2:end);
    
    dline = spdiags(dxdxiNorm.*IPw_f,0,g,g);
    M = (referenceElement.N1d)'*(dline*referenceElement.N1d);
    b1 = (referenceElement.N1d)'*(dline*ux_hatI);
    b2 = (referenceElement.N1d)'*(dline*uy_hatI);
    uf_x = M\b1;
    uf_y = M\b2;
    uface = [uf_x ; uf_y];

    %press_I_e=analyticalpressure(Pgauss,time);
    %
    if(cnt==1)
        press_I_e = press_initial_condition(Pgauss,time);
    else
        press_I_e = analyticalpressure(Pgauss,time);
    end 
    %
    dlineP = spdiags(dxdxiNorm.*IPw_f,0,g,g);
    M_P = (referenceElement.N1d)'*(dlineP*referenceElement.N1d);
    b_P = (referenceElement.N1d)'*(dlineP*press_I_e);    
    press_face = M_P\b_P;
    
    %u_I = [uf_x , uf_y];
      
         %figure(1); hold on;
         %plot(Pphy(:,1),Pphy(:,2),'ro');
         %quiver(I(:,1),I(:,2),Iprime(:,1),Iprime(:,2));
    
    fu_I = (AuuhatI*uface);
    fp_I = (AuphatI*press_face);
    fq_I = (AquhatI*uface);
    
    % Elemental mapping for terms containing surface integrals
    Alq = -Aql'; Alu = Aul';
    %%% compute Auu from it's different components
    Auu = Auu_vol_t+Auu_surf+AuuI;

    %%% compute fe from it's different components
    fe = fe_s_vel + fe_uns - fe_s_press_surf + fe_s_press_vol + fu_I - fp_I;
    %fe_2 = fe_s_vel + fe_uns - fe_s_press_surf + fe_s_press_vol + (AuuhatI*uhatI_e_old);
    fq = fq_I;
    %fq_2 = AquhatI*uhatI_e_old;
    %%% compute the source term for the global eqn.
    %g = gn_press;
    %
    %Alq = -Aql'; Alu = Aul';
    A = [Auu Auq; Aqu Aqq];
    UQ = -A\[Aul;Aql];
    %fUQ= A\[(fe-(AuuhatI*uhatI_e));-(AquhatI*uhatI_e)];
    fUQ= A\[fe;fq];
    %U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
    %Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q
    U = UQ(1:2*nOfElementNodes,:); Uf=fUQ(1:2*nOfElementNodes); % maps lamba into U
    Q = UQ(2*nOfElementNodes+1:end,:); Qf=fUQ(2*nOfElementNodes+1:end); % maps lamba into Q
    disp('Hola');
        
   %{
     %%  *****LOCAL TEST FOR CUT ELEMENT*****
     %u_analy=analiticalSolutionLaplace(Xe,time);
     u_analy = analyticalVelocityStokes(Xe,time);
     x = Xe(:,1);
     y = Xe(:,2);
     %
     %
     %dux_dx = 0+0.*x; duy_dx = 0+0.*y; dux_dy = 0+0.*x; duy_dy = 0+0.*y;
     %dux_dx = 1+0.*x; duy_dx = 0+0.*y; dux_dy = 0+0.*x; duy_dy = 0+0.*y;
     dux_dx = 1+0.*x; duy_dx = 1+0.*y; dux_dy = 1+0.*x; duy_dy = -1+0.*y;    
     %dux_dx = 0+2.*x; dux_dy = 0+0.*x; duy_dx = 0-2.*y; duy_dy = 0-2.*x;
     %dux_dx = 2.*x; dux_dy = 0.*x; duy_dx = -1.*y; duy_dy = -1.*x;
     %dux_dx = -2.*y; dux_dy = -2.*x; duy_dx = 0.*y; duy_dy = 2.*y;
     %dux_dx = -1.*y; dux_dy = -1.*x; duy_dx = 0.*y; duy_dy = 2.*y;
     %dux_dx = 2.*x.*y; dux_dy = x.^2; duy_dx = -1.*y.^2; duy_dy = -2.*x.*y;
     %dux_dx = 3.*x.*x; dux_dy = 0.*x; duy_dx = -6.*x.*y; duy_dy = -3.*x.*x;
     %dux_dx = -3.*y.*y; dux_dy = -6.*x.*y; duy_dx = 0.*x; duy_dy = 3.*y.*y;
     %dux_dx = time+0.*x; duy_dx = 0+0.*y; dux_dy = 0+0.*x; duy_dy = time+0.*y;
     %dux_dx = (-2*time).*y; duy_dx = (-2*time).*x; dux_dy = 0+0.*x; duy_dy = (2*time).*y;
     %
     q_xx = dux_dx./Re; q_yx = duy_dx./Re; q_xy = dux_dy./Re; q_yy = duy_dy./Re;
     %
     q_analy = [q_xx; q_xy; q_yx; q_yy];
%     
aux2=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1); 
gama=[analyticalVelocityStokes(Xe(aux2(1:nOfFaceNodes),:),time);analyticalVelocityStokes(Xe(aux2(nOfFaceNodes+1:2*nOfFaceNodes),:),time);...
   analyticalVelocityStokes(Xe(aux2(2*nOfFaceNodes+1:3*nOfFaceNodes),:),time)];

%{         
     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama_x=Xe(ce,1);
     coord_gama_y=Xe(ce,2);
     coord_gama=[coord_gama_x,coord_gama_y];
     gama=analiticalSolutionLaplace(coord_gama,time);
%}     
     uhat_analy= analyticalVelocityStokes(Pphy,time);
     phat_analy= analyticalpressure(Pphy,time);
     a1 = Aqu*u_analy;
     a2 = Aqq*q_analy;
     a3 = Aql*gama;
     a4 = AquhatI*uhat_analy;
     tetstest=Aqu*u_analy+Aqq*q_analy+Aql*gama-AquhatI*uhat_analy;
%     
     disp(['Cut Element: TesttestElem' num2str(iElem) ]);
     disp(max(abs(tetstest)));
   
    
    %a=AquhatI*uhat_analy;
    %disp(max(a));
    
% %     
% %     %smalltest
% %
     %{
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
     %}
% %    
     %
     a_1 = Auu*u_analy;
     b_1 = Aul*gama;
     c_1 = Auq*q_analy;
     d_1 = AuuhatI*uhat_analy;
     %
     fe_temp = fe_s_vel+fe_uns-fe_s_press_surf+fe_s_press_vol;
% %     
     test_result_1=fe_temp-(Auu*u_analy)-(Auq*q_analy)-(Aul*gama)+(AuuhatI*uhat_analy)-(AuphatI*phat_analy);
     test_result_2=(Aqq*q_analy)+(Aqu*u_analy)+(Aql*gama)-(AquhatI*uhat_analy);
% %     
     disp(['Cut element: Test result equation 1 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_1)));
     disp(['Cut element: Test result equation 2 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_2)));
     disp('Hola')
     %}
% %     
% % 
%
    
elseif ~isempty(c) % for standard HDG element
    
    % Elemental mapping
    %%% compute fe from it's different components
    fe = fe_s_vel + fe_uns - fe_s_press_surf + fe_s_press_vol;
    %
    fq = zeros(4*nOfElementNodes,1);
    %%% compute the source term for the global eqn.
    %g = gn_press;
    %
    Auu = Auu_vol_t+Auu_surf; 
    Alq = -Aql'; Alu = Aul';
    A = [Auu Auq; Aqu Aqq];
    UQ = -A\[Aul;Aql];
    %fUQ= A\[(fe_s+fe_uns);zeros(2*nOfElementNodes,1)];
    fUQ = A\[fe;fq];
    %U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
    %Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q
    U = UQ(1:2*nOfElementNodes,:); Uf=fUQ(1:2*nOfElementNodes); % maps lamba into U
    Q = UQ(2*nOfElementNodes+1:end,:); Qf=fUQ(2*nOfElementNodes+1:end); % maps lamba into Q     

    %{
% %     %small test
%% ********************************************************************
%%      %LOCAL PROBLEM FOR STANDARD HDG ELEMENT 
% % %     % % Analytical Solution     
% %
%
%
%
aux2=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1); 
gama=[analyticalVelocityStokes(Xe(aux2(1:nOfFaceNodes),:),time);analyticalVelocityStokes(Xe(aux2(nOfFaceNodes+1:2*nOfFaceNodes),:),time);...
   analyticalVelocityStokes(Xe(aux2(2*nOfFaceNodes+1:3*nOfFaceNodes),:),time)];


     %%% u_analy = analiticalSolutionLaplace(Xe);
     u_analy = analyticalVelocityStokes(Xe,time);
     %%% q_analy=zeros(2*length(u_analy),1); q_analy(1:2:2*nOfElementNodes)=-25*Xe(:,1).^4;
     %aux = zeros(nOfElementNodes,1);
     x = Xe(:,1);
     y = Xe(:,2);
     %
     %
     %dux_dx = 0+0.*x; duy_dx = 0+0.*y; dux_dy = 0+0.*x; duy_dy = 0+0.*y;
     dux_dx = 1+0.*x; duy_dx = 1+0.*y; dux_dy = 1+0.*x; duy_dy = -1+0.*y;    
     %dux_dx = 0+2.*x; dux_dy = 0+0.*x; duy_dx = 0-2.*y; duy_dy = 0-2.*x;
     %dux_dx = 2.*x; dux_dy = 0.*x; duy_dx = -1.*y; duy_dy = -1.*x;
     %dux_dx = -2.*y; dux_dy = -2.*x; duy_dx = 0.*y; duy_dy = 2.*y;
     %dux_dx = -1.*y; dux_dy = -1.*x; duy_dx = 0.*y; duy_dy = 2.*y;
     %dux_dx = 2.*x.*y; dux_dy = x.^2; duy_dx = -1.*y.^2; duy_dy = -2.*x.*y;
     %dux_dx = 3.*x.*x; dux_dy = 0.*x; duy_dx = -6.*x.*y; duy_dy = -3.*x.*x;
     %dux_dx = -3.*y.*y; dux_dy = -6.*x.*y; duy_dx = 0.*x; duy_dy = 3.*y.*y;
     %dux_dx = time+0.*x; duy_dx = 0+0.*y; dux_dy = 0+0.*x; duy_dy = time+0.*y;
     %dux_dx = (-2*time).*y; duy_dx = (-2*time).*x; dux_dy = 0+0.*x; duy_dy = (2*time).*y;
     %
     q_xx = dux_dx./Re; q_yx = duy_dx./Re; q_xy = dux_dy./Re; q_yy = duy_dy./Re;
     %
     q_analy = [q_xx; q_xy; q_yx; q_yy];
     %
     %
     p_analy = 0.*(x+y);
     %
     %
% %     
     %%% test_result_1=fe-(Auu*u_analy)-(Auq*q_analy)-(Aul*gama);
     %%% test_result_2=(Aqq*q_analy)+(Aqu*u_analy)+(Aql*gama);
     %
     a_1 = Auu*u_analy;
     b_1 = Aul*gama;
     c_1 = Auq*q_analy;
     %
     %
     test_result_1 = fe-(Auu*u_analy)-(Auq*q_analy)-(Aul*gama);
     %
     %
     aa = Aqu*u_analy;
     bb = Aqq*q_analy;
     cc = Aql*gama;
% % % %      %dd = Aul*gama;
% % % % 
     test_result_2 = (Aqu*u_analy)+(Aqq*q_analy)+(Aql*gama);
% %     
     disp(['Standard (UnCut) Element: Test result equation 1 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_1)));
     disp(['Standard (UnCut) Element: Test result equation 2 for element ' num2str(iElem) ]);
     disp(max(abs(test_result_2)));
     disp('Hola')
     %}
  
elseif ~isempty(d) % for elements in void region
    
    U =zeros(2*nOfElementNodes,3*2*nOfFaceNodes);
    Q =zeros(4*nOfElementNodes,3*2*nOfFaceNodes);
    Uf=zeros(2*nOfElementNodes,1);
    Qf=zeros(4*nOfElementNodes,1);
    
end

%disp('Hola')













