function [KK,f, QQ, UU, Qf, Uf, uDirichlet] = hdg_matrix_predictorStep(X,T,F,referenceElement,infoFaces,tau,dt,cnt,u_old,u_old2,u_old3,p_old,time,Re,a_parm)
% This routine is valid only for triangles

nOfElements = size(T,1);
nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = max(max(F));
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

aux=nOfFaceNodes:-1:1; aux=[aux,aux+nOfFaceNodes]; 
indflip=[aux,2*nOfFaceNodes+aux,4*nOfFaceNodes+aux];
%u_old = ProjectingOldSolution(u,nOfElements,X,T,referenceElement);

% loop in elements
for iElem =1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)
    %mu=mu_vector(iElem);
    %a_cnt = (iElem-1)*2*nOfElementNodes+(1:2*nOfElementNodes);
    a_cnt = (iElem-1)*nOfElementNodes+(1:nOfElementNodes);
    %
    u_intm1 = u_old(a_cnt,:);
    u_Elem_old = [u_intm1(:,1) ; u_intm1(:,2)];
    %
    u_intm2 = u_old2(a_cnt,:);
    u_Elem_old2 = [u_intm2(:,1) ; u_intm2(:,2)];
    %
    u_intm3 = u_old3(a_cnt,:);
    u_Elem_old3 = [u_intm3(:,1) ; u_intm3(:,2)];

    %
    p_Elem = p_old(a_cnt);
    %u_Elem = [u(a_cnt,1) ; u(a_cnt,2)];

    % elemental matrices
    [Qe,Ue,Qfe,Ufe,Alq,Alu,All,g] = KKeElementalMatrices(Xe,referenceElement,tau,dt,cnt,u_Elem_old,u_Elem_old2,u_Elem_old3,p_Elem,time,Re,a_parm);
    
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
    ffe = g-(Alq*Qfe + Alu*Ufe);
    
    %All faces are assembled (Dirichlet included)
    aux = (1:2*nOfFaceNodes);
    ind = [(Fe(1)-1)*2*nOfFaceNodes + aux,(Fe(2)-1)*2*nOfFaceNodes + aux,(Fe(3)-1)*2*nOfFaceNodes + aux];
    KK(ind,ind)=KK(ind,ind)+ KKe;
    f(ind) = f(ind) + ffe;
    %
    %disp('Hola')
    %
end

%Dirichlet BC
%Dirichlet face nodal coordinates
nOfExteriorFaces=size(infoFaces.extFaces,1);
XDirichlet = zeros(nOfExteriorFaces*nOfFaceNodes,2);
for iFace=1:nOfExteriorFaces
    iElem = infoFaces.extFaces(iFace,1);
    nFace = infoFaces.extFaces(iFace,2);
    ind = (iFace-1)*nOfFaceNodes+[1:nOfFaceNodes];
    XDirichlet(ind,:) = X(T(iElem,referenceElement.faceNodes(nFace,:)),:);
end
figure(1), hold on, plot(XDirichlet(:,1),XDirichlet(:,2),'r*'), hold off
%uDirichlet = DirichletCondition(XDirichlet);
%uDirichlet = DirichletCondition(XDirichlet);
%udbc = @DirichletCondition;
uDirichlet=computeProjectionFacesDirchletPredictorStep(infoFaces.extFaces,X,T,referenceElement,time,Re);
% System reduction (Dirichlet faces are set to prescribed value)
%nDOF = nOfInteriorFaces*nOfFaceNodes;
nDOF = 2*nOfInteriorFaces*nOfFaceNodes;
f = f(1:nDOF)-KK(1:nDOF,nDOF+1:end)*uDirichlet;
KK=KK(1:nDOF,1:nDOF);
%
%disp('Hola')
%

%%
%% ELEMENTAL MATRIX
function [Q,U,Qf,Uf,Alq,Alu,All,gn_press] = KKeElementalMatrices(Xe,referenceElement,tau,dt,cnt,u_Elem_old,u_Elem_old2,u_Elem_old3,p_Elem,time,Re,a_parm)
                                                        
nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)
%tau = [tau,zeros(1,nOfFaces-1)]; %Non-null only for 1st face in each element

% Information of the reference element
N = referenceElement.N;
Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;

%Numerical quadrature
IPw = referenceElement.IPweights; ngauss=size(IPw,1);
%Xg = N*Xe; %physical coordinates
IPw_f = referenceElement.IPweights1d;

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);
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

% call the subroutine which computes the source term.
%sourceTerm = @analiticalLaplacianLaplace;

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
    dvolu=IPw(g)*det(J);
    
    %physical coordinates of the Gauss points on current Element
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
    
    %{
    if(cnt==1)
        u_old_gauss = vel_initial_condition(xy_g);
    else
        u_old_gauss = cg_2*u_Elem; 
    end 
    %
    %fe_uns = fe_uns + N_g'*u_old_gauss*dvolu/dt;
    fe_uns = fe_uns + cg_1*u_old_gauss*dvolu/(a_parm*dt);
    %}
    %
    if(cnt==1)        %for n=1
        u_old_gauss = vel_initial_condition(xy_g);
        %here: u_old_gauss = u_(n-1)
    elseif(cnt==2)    %for n=2
        u_old_gauss = vel_initial_condition(xy_g);
        u_old_gauss2 = cg_2*u_Elem_old2;
        %here: u_old_gauss2 = u_(n-1) and u_old_gauss = u_(n-2)
    elseif(cnt==3)    %for n=3
        u_old_gauss = vel_initial_condition(xy_g);
        u_old_gauss2 = cg_2*u_Elem_old2;       
        u_old_gauss3 = cg_2*u_Elem_old3;
        %here: u_old_gauss3 = u_(n-1) and u_old_gauss2 = u_(n-2)
        % and u_old_gauss = u_(n-3)
    else              %for n>3
        u_old_gauss = cg_2*u_Elem_old;
        u_old_gauss2 = cg_2*u_Elem_old2;
        u_old_gauss3 = cg_2*u_Elem_old3;
        %here: u_old_gauss3 = u_(n-1) and u_old_gauss2 = u_(n-2)
        % and u_old_gauss = u_(n-3)
    end 
    %u_old_gauss = N_g*u_Elem;
    
    if(cnt==1)
        %
        fe_uns = fe_uns + cg_1*u_old_gauss*dvolu/dt;
    elseif(cnt==2)
        %
        aa = -1*cg_1*u_old_gauss*dvolu/(2*dt);
        bb = 2*cg_1*u_old_gauss2*dvolu/dt;
        fe_uns = fe_uns + aa + bb;
    elseif(cnt==3)
        %
        aa = cg_1*u_old_gauss*dvolu/(3*dt);
        bb = -3*cg_1*u_old_gauss2*dvolu/(2*dt);
        cc = 3*cg_1*u_old_gauss3*dvolu/dt;
        fe_uns = fe_uns + aa + bb + cc;
    else
        %
        aa = cg_1*u_old_gauss*dvolu/(3*dt);
        bb = -3*cg_1*u_old_gauss2*dvolu/(2*dt);
        cc = 3*cg_1*u_old_gauss3*dvolu/dt;
        fe_uns = fe_uns + aa + bb + cc;
    end
    
    
end
%inv_Aqq = inv(Aqq);
% Elemental mapping for terms containing volume integrals
Aqu = -Auq';

%%%***Allocate the arrays Aii used in surf. integrals of local elemental problem***
Auu_surf = zeros(2*nOfElementNodes,2*nOfElementNodes);
Alq = zeros(2*3*nOfFaceNodes,4*nOfElementNodes);
Aql = zeros(4*nOfElementNodes,2*3*nOfFaceNodes);
Alu = zeros(2*3*nOfFaceNodes,4*nOfElementNodes);
Aul = zeros(2*nOfElementNodes,2*3*nOfFaceNodes);
All = zeros(2*3*nOfFaceNodes,2*3*nOfFaceNodes);
gn_press = zeros(2*3*nOfFaceNodes,1);

%% Faces computations:
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces
    
    tau_f = tau(iface);
    
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
    
    % Gauss points position on the current face
    xfg = N1d*xf;
    yfg = N1d*yf;   
  
    % Inizialization 
    %allocate local arrays for local computations on the faces.
    Auu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    Aql_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    Aql_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_f = zeros(nOfFaceNodes,nOfFaceNodes);
    All_f = zeros(nOfFaceNodes,nOfFaceNodes);   
    %Alu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    fe_s_press_fx = zeros(nOfFaceNodes,1);
    fe_s_press_fy = zeros(nOfFaceNodes,1);
      
    %  LOOP IN GAUSS POINTS
    for g = 1:ngauss_f
        
        % Shape functions and derivatives at the current integration point
        Nf_g = N1d(g,:);
        Nfxi_g = Nx1d(g,:);
        
        % Integration weight
        xyDer_g = Nfxi_g*[xf yf];
        xyDerNorm_g = norm(xyDer_g);       
        dline=IPw_f(g)*xyDerNorm_g;
        pp = sqrt(xyDer_g(1).^2+xyDer_g(2).^2);
        
        % Unit normal to the boundary
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) -t_g(1)];
        
        % physical location of the Gauss points
        xy_gf = [xfg(g) yfg(g)];
        %
        % Contribution of the current integration point to the elemental matrix
        %
         Auu_f = Auu_f + tau_f*(Nf_g')*Nf_g*dline;
         Aql_fnx = Aql_fnx - Nf_g'*Nf_g*n_g(1)*dline;
         Aql_fny = Aql_fny - Nf_g'*Nf_g*n_g(2)*dline;
         Aul_f = Aul_f - tau_f*(Nf_g')*Nf_g*dline;
         All_f = All_f + tau_f*(Nf_g')*Nf_g*dline;
         
         %compute the second pressure source term 
         if(cnt==1)
             p_old_gaussFace = press_initial_condition(xy_gf,time);
         else
             p_old_gaussFace = Nf_g*p_Elem(nodes);
         end
         dd_x = (Nf_g'*n_g(1))*p_old_gaussFace*dline;
         dd_y = (Nf_g'*n_g(2))*p_old_gaussFace*dline;
         % for the use in local eqn. of predictor velocity (in source term)
         % and for the use in the global eqn. (in source term)
         fe_s_press_fx = fe_s_press_fx + dd_x;
         fe_s_press_fy = fe_s_press_fy + dd_y; 
         %
    end
    %    
    Auu_surf(nodes1,nodes1) = Auu_surf(nodes1,nodes1)+Auu_f; 
    Auu_surf(nodes2,nodes2) = Auu_surf(nodes2,nodes2)+Auu_f; 

    Aul(nodes1,ind_face1) = Aul(nodes1,ind_face1)+Aul_f; 
    Aul(nodes2,ind_face2) = Aul(nodes2,ind_face2)+Aul_f;

    Aql(nodes1,ind_face1) = Aql(nodes1,ind_face1)+Aql_fnx;
    Aql(nodes2,ind_face1) = Aql(nodes2,ind_face1)+Aql_fny;
    Aql(nodes3,ind_face2) = Aql(nodes3,ind_face2)+Aql_fnx;
    Aql(nodes4,ind_face2) = Aql(nodes4,ind_face2)+Aql_fny; 
    
    All(ind_face1,ind_face1)=All(ind_face1,ind_face1)+All_f;
    All(ind_face2,ind_face2)=All(ind_face2,ind_face2)+All_f;
    
    fe_s_press_surf(nodes1) = fe_s_press_surf(nodes1)+fe_s_press_fx;
    fe_s_press_surf(nodes2) = fe_s_press_surf(nodes2)+fe_s_press_fy;
    
    gn_press(ind_face1) = gn_press(ind_face1)+fe_s_press_fx;
    gn_press(ind_face2) = gn_press(ind_face2)+fe_s_press_fy;
    
end

% Elemental mapping for terms containing surface integrals
Alq = -Aql'; Alu = Aul';

%%% compute Auu from it's different components
Auu = Auu_vol_t + Auu_surf; 
%%% compute fe from it's different components
fe = fe_s_vel + fe_uns - fe_s_press_surf + fe_s_press_vol;
%%% compute the source term for the global eqn.
%g = gn_press;

A = [Auu Auq; Aqu Aqq];
UQ = -A\[Aul;Aql];

fq = zeros(4*nOfElementNodes,1);
fUQ= A\[fe;fq];

U = UQ(1:2*nOfElementNodes,:); Uf=fUQ(1:2*nOfElementNodes); % maps lamba into U
Q = UQ(2*nOfElementNodes+1:end,:); Qf=fUQ(2*nOfElementNodes+1:end); % maps lamba into Q 

%disp('Say Hola')
%% ********************************************************************
% % %     %TEST OF LOCAL PROBLEM FOR THE PREDICTOR VELOCITY 
% % %     % % Analytical Solution     
% %
%
%{
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
     %
     %
     if any([max(abs(test_result_1)),max(abs(test_result_2))]>1.e-14)
             disp(sprintf('      LP max(abs(residual)) eq. #1: %e',max(abs(test_result_1))));
             disp(sprintf('      LP max(abs(residual)) eq. #2: %e',max(abs(test_result_2))));
     else
         disp('LP OK')
     end

     disp('Hola')
%}
%
%
%% **************
%{
% % %     %TEST OF LOCAL PROBLEM
% % %     % % Analytical Solution     
     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama=Xe(ce,:);
% % 
     u_analy = analiticalSolutionLaplace(Xe,time);
     q_analy=zeros(2*length(u_analy),1); 
     
     x = Xe(:,1);
     y = Xe(:,2);

     %trmA = sin(pi*x)+pi.*cos(pi*x);
     du_dx = 1+0.*x;
     %q_analy(1:2:2*nOfElementNodes)=-25*Xe(:,1).^4;
     q_analy(1:2:2*nOfElementNodes)= -1.*du_dx;

     %trmB = sin(pi*y)+pi.*cos(pi*y);
     du_dy = 0+0.*y;
     q_analy(2:2:2*nOfElementNodes)= -1.*du_dy;
     %q_analy = [-1.*du_dx;-1.*du_dy];

     gama = analiticalSolutionLaplace(coord_gama,time);

     %a1 = Auu_vol*u_analy;
     %a2 = Auu_surf*u_analy;
     %b1 = Aul_1*gama;
     a = Auu*u_analy;
     b = Aul*gama;
     c = Auq*q_analy;
% %     
     test_result_1=fe_s+fe_uns-(Auu*u_analy)-(Auq*q_analy)-(Aul*gama);
     test_result_2=(Aqq*q_analy)+(Aqu*u_analy)+(Aql*gama);
% %     
     if any([max(abs(test_result_1)),max(abs(test_result_2))]>1.e-14)
             disp(sprintf('      LP max(abs(residual)) eq. #1: %e',max(abs(test_result_1))));
             disp(sprintf('      LP max(abs(residual)) eq. #2: %e',max(abs(test_result_2))));
     else
         disp('LP OK')
     end

     disp('Hola')
%}


% %     u_analy=analiticalSolutionLaplace(Xe);
% %     x = Xe(:,1);
% %     y = Xe(:,2);
% %     
% % 
% %     
% %     q_analy=zeros(2*length(u_analy),1);
% %     q_analy(1:2:end)=-1*length(u_analy);
% %     q_analy(2:2:end)=0*length(u_analy);
% %     
% %     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
% %     coord_gama_x=Xe(ce,1);
% %     coord_gama_y=Xe(ce,2);
% %     coord_gama=[coord_gama_x,coord_gama_y];
% %     gama=analiticalSolutionLaplace(coord_gama);
% %     
% %     
% %     
% %     test_result_1=fe-(Auu*u_analy)-(Auq*q_analy)-(Aul*gama);
% %     test_result_2=(Aqq*q_analy)+(Aqu*u_analy)+(Aql*gama);
% %     
% %     
% %     if abs(test_result_1)>10^-10
% %     disp(['Test result equation 1 for element ' num2str(iElem) ]);
% %     disp(max(abs(test_result_1)));
% %     end
% %     
% %     if abs(test_result_2)>10^-10
% %     disp(['Test result equation 1 for element ' num2str(iElem) ]);
% %     disp(max(abs(test_result_2)));
% %     end
% %         
    
% %     disp(['Test result equation 1 for element ' num2str(iElem) ]);
% %     disp(max(abs(test_result_1)));
% %     disp(['Test result equation 2 for element ' num2str(iElem) ]);
% %     disp(max(abs(test_result_2)));
% % 
    %Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Aql_fnx;
    %Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Aql_fny;
    %Alu(ind_face,nodes) = Alu(ind_face,nodes) + Alu_f;
    %Alu(ind_face,nodes) = Alu(ind_face,nodes) + Auu_f;
    
    %Aql(2*nodes-1,ind_face) = Aql(2*nodes-1,ind_face) + Aql_fnx;
    %Aql(2*nodes,ind_face) = Aql(2*nodes,ind_face) + Aql_fny;
    %Aul(nodes,ind_face) = Aul(nodes,ind_face) + Aul_f;    
    %Auu_surf(nodes,nodes) = Auu_surf(nodes,nodes) + Auu_f;    
    %All(ind_face,ind_face) = All(ind_face,ind_face) + All_f;
         
                    
        % Auu_f = Auu_f + tau_f *mu*(Nf_g')*Nf_g*dline;
        %%% All_f = All_f - tau_f *2*mu*(Nf_g')*Nf_g*dline;
        %All_f = All_f - tau_f *2*(Nf_g')*Nf_g*dline;
        %Alq_fnx = Alq_fnx + Nf_g'*Nf_g*n_g(1)*dline;
        %Alq_fny = Alq_fny + Nf_g'*Nf_g*n_g(2)*dline;
        %Alu_f = Alu_f + tau_f*(Nf_g')*Nf_g*dline; 
        %Aul_f = Aul_f + bb*(Nf_g)'*Nf_g*dline;

%Auu_vol_t = zeros(nOfElementNodes,nOfElementNodes);
%Auu_vol_s = zeros(nOfElementNodes,nOfElementNodes);
%Auq = zeros(nOfElementNodes,2*nOfElementNodes);
%Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
%fe_s = zeros(1*nOfElementNodes,1);
%fe_s_press = zeros(2*nOfElementNodes,1);
%fe_uns = zeros(1*nOfElementNodes,1);



