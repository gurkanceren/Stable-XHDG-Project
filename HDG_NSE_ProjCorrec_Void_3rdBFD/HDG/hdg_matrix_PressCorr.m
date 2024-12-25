function [KK,f, QQ, UU, Qf, Uf, pDirichlet] = hdg_matrix_PressCorr(X,T,F,referenceElement,infoFaces,tauP,dt,Re,a_parm,u,uhat,time)
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
for iElem =1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
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
    [Qe,Ue,Qfe,Ufe,Alq,Alp,All] = KKeElementalMatrices(Fe,Xe,referenceElement,nOfInteriorFaces,infoFaces,u_Elem,uhat_Elem,tauP,dt,a_parm,iElem,time);
    
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
    %KKe = Alq*Qe + Alp*Ue + 0.5*All;
    KKe = Alq*Qe + Alp*Ue + All;
    ffe = -(Alq*Qfe + Alp*Ufe);
    
    %All faces are assembled (Dirichlet included)
    aux = (1:nOfFaceNodes);
    ind = [(Fe(1)-1)*nOfFaceNodes + aux,(Fe(2)-1)*nOfFaceNodes + aux,(Fe(3)-1)*nOfFaceNodes + aux];
    KK(ind,ind)=KK(ind,ind)+ KKe;
    f(ind) = f(ind) + ffe;
    
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
%pDirichlet=computeProjectionFacesDirchlet(@DirichletCondition,infoFaces.extFaces,X,T,referenceElement);
pDirichlet=computeProjectionFacesDirchletPressCorr(infoFaces.extFaces,X,T,referenceElement,time,Re);
% System reduction (Dirichlet faces are set to prescribed value)
nDOF = nOfInteriorFaces*nOfFaceNodes;
f = f(1:nDOF)-KK(1:nDOF,nDOF+1:end)*pDirichlet;
KK=KK(1:nDOF,1:nDOF);
%
%disp('Hola_Hola')

%%
%% ELEMENTAL MATRIX
function [Q,U,Qf,Uf,Alq,Alp,All] = KKeElementalMatrices(Fe,Xe,referenceElement,nOfInteriorFaces,infoFaces,u_Elem,uhat_Elem,tau,dt,a_parm,iElem,time)

nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)
%tau = [tau,zeros(1,nOfFaces-1)]; %Non-null only for 1st face in each element

uhat_Elem_X = zeros(nOfElementNodes,1);
uhat_Elem_Y = zeros(nOfElementNodes,1);
fe = zeros(nOfElementNodes,1);
fe_vol = zeros(nOfElementNodes,1);
fe_surf_1 = zeros(nOfElementNodes,1);
fe_surf_2 = zeros(nOfElementNodes,1);
Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
Apq = zeros(nOfElementNodes,2*nOfElementNodes);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
App = zeros(nOfElementNodes,nOfElementNodes);
Alp = zeros(3*nOfFaceNodes,nOfElementNodes);
Apl = zeros(nOfElementNodes,3*nOfFaceNodes);
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

%%% call the subroutine which computes the source term.
%%% sourceTerm = @analiticalLaplacianLaplace;

%% Volume computation: loop in Gauss points
for g = 1:ngauss
    
    %Shape functions and derivatives at the current Gauss point
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
    
    %physical coordinates of the Gauss points on current Element
    xy_g=N_g*Xe;
    
    %x and y derivatives
    invJ = inv(J);  
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    NN_xy(1:2:end) = Nx_g;    NN_xy(2:2:end) = Ny_g;
    
    %Contribution of the current Gauss point to the elemental matrix
    %
    Aqq = Aqq + NN_g'*NN_g*dvolu;
    Apq = Apq - N_g'*NN_xy*dvolu;
    %
    grad_testfunc = [Nx_g ; Ny_g];
    grad_testfunc_T = grad_testfunc';
    %
    u_gauss = N_g*u_Elem;
    %
    fe_vol = fe_vol +grad_testfunc_T*u_gauss'*dvolu/(a_parm*dt); %compute first source term
    %
end
%inv_Aqq = inv(Aqq);

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
    
    tau_f = tau(iface);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    %
    % local numbering of node points on mesh skeleton faces.
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    %ind_face1 = (iface-1)*2*nOfFaceNodes + (1:nOfFaceNodes); %% for x comp.
    %ind_face2 = ind_face1 + nOfFaceNodes;  %% for y comp.
    % 
    xf = xe(nodes);
    yf = ye(nodes);
    
    % Gauss points position on the current face
    xfg = N1d*xf;
    yfg = N1d*yf;   
  
    % Inizialization
    Alq_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    Alq_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    Alp_f = zeros(nOfFaceNodes,nOfFaceNodes);
    Apl_f = zeros(nOfFaceNodes,nOfFaceNodes);
    App_f = zeros(nOfFaceNodes,nOfFaceNodes);
    All_f = zeros(nOfFaceNodes,nOfFaceNodes);
    fe_uhat_fx = zeros(nOfFaceNodes,1);
    fe_uhat_fy = zeros(nOfFaceNodes,1);
    fe_u_fx = zeros(nOfFaceNodes,1);
    fe_u_fy = zeros(nOfFaceNodes,1);
    
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
        %
        %
        % Contribution of the current Gauss point to the elemental matrices
        %   
        %All_f = All_f - tau_f *2*(Nf_g')*Nf_g*dline;
        All_f = All_f + tau_f*(Nf_g')*Nf_g*dline;
        %
        Alq_fnx = Alq_fnx + Nf_g'*Nf_g*n_g(1)*dline;
        Alq_fny = Alq_fny + Nf_g'*Nf_g*n_g(2)*dline;
        %
        %Alu_f = Alu_f + tau_f*(Nf_g')*Nf_g*dline; 
        %
        App_f = App_f + tau_f*(Nf_g')*Nf_g*dline;
        %
        Apl_f = Apl_f - tau_f*(Nf_g)'*Nf_g*dline;
        %
        %compute the second source term
        %uhatX_gaussFace = Nf_g*uhat_Elem(ind_face1);
        %uhatY_gaussFace = Nf_g*uhat_Elem(ind_face2); 
        %
        uhatX_gaussFace = Nf_g*uhat_Elem_XN(ind_face);
        uhatY_gaussFace = Nf_g*uhat_Elem_YN(ind_face);        
        %
        dd_hat_x = (Nf_g'*n_g(1))*uhatX_gaussFace*dline/(a_parm*dt);
        dd_hat_y = (Nf_g'*n_g(2))*uhatY_gaussFace*dline/(a_parm*dt);
        %
        %
        fe_uhat_fx = fe_uhat_fx - dd_hat_x;
        fe_uhat_fy = fe_uhat_fy - dd_hat_y;
        %
        %
        uX_gaussFace = Nf_g*u_Elem(nodes,1);
        uY_gaussFace = Nf_g*u_Elem(nodes,2);
        %
        dd_x = (Nf_g'*n_g(1))*uX_gaussFace*dline/(a_parm*dt);
        dd_y = (Nf_g'*n_g(2))*uY_gaussFace*dline/(a_parm*dt);
        %
        % 
        fe_u_fx = fe_u_fx + dd_x;
        fe_u_fy = fe_u_fy + dd_y;
        %          
    end
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


% Elemental mapping
Aqp = -Apq'; Aql = -Alq'; Alp = Apl'; 
%
fe = fe_vol + fe_surf_1 +fe_surf_2;


A = [App Apq; Aqp Aqq];
UQ = -A\[Apl;Aql];

fq = zeros(2*nOfElementNodes,1);
fUQ= A\[fe;fq];

U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q

%disp('Hola')
%
%{
% % %     %TEST OF LOCAL PROBLEM
% % %     % % Analytical Solution     
     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama=Xe(ce,:);
% % 
     
     %u_analy = analiticalSolutionLaplace(Xe);
     p_analy = analyticalpressure(Xe,time);
     q_analy=zeros(2*length(p_analy),1); 
     
     x = Xe(:,1);
     y = Xe(:,2);

     %trmA = sin(pi*x)+pi.*cos(pi*x);
     dp_dx = time+0.*x;
     %q_analy(1:2:2*nOfElementNodes)=-25*Xe(:,1).^4;
     %q_analy(1:2:2*nOfElementNodes)= -1.*du_dx;
     q_analy(1:2:2*nOfElementNodes)= dp_dx;

     %trmB = sin(pi*y)+pi.*cos(pi*y);
     dp_dy = time+0.*y;
     %q_analy(2:2:2*nOfElementNodes)= -1.*du_dy;
     q_analy(2:2:2*nOfElementNodes)= dp_dy;
     %q_analy = [-1.*du_dx;-1.*du_dy];

     %gama = analiticalSolutionLaplace(coord_gama);
     gama = analyticalpressure(coord_gama,time);

     %
     %a1 = App_vol*u_analy;
     %a2 = App*u_analy;
     %b1 = Apl_1*gama;
     aa_1 = App*p_analy;
     aa_2 = Apl*gama;
     aa_3 = Apq*q_analy;
     %
% %     
     test_result_1=fe-(App*p_analy)-(Apq*q_analy)-(Apl*gama);
     test_result_2=(Aqq*q_analy)+(Aqp*p_analy)+(Aql*gama);
% %     
     if any([max(abs(test_result_1)),max(abs(test_result_2))]>1.e-14)
             disp(sprintf('      LP max(abs(residual)) eq. #1: %e',max(abs(test_result_1))));
             disp(sprintf('      LP max(abs(residual)) eq. #2: %e',max(abs(test_result_2))));
     else
         disp('LP OK')
     end

     disp('Lola')
%}
%

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






