function [u_crc,p_crc,p_wrc_crc]=computeElementalCorrectionStep_Qcorr(u,uhat,pcorr,phat,qcorr,X,T,F,referenceElement,infoFaces,tauP,dt,Re,a_parm,time,cnt,pold)
% This routine is valid only for triangles
%global filter;

nOfElements = size(T,1);
nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = max(max(F));
%
u_crc = zeros(nOfElements*nOfElementNodes,2);
p_crc = zeros(nOfElements*nOfElementNodes,1);
p_wrc_crc = zeros(nOfElements*nOfElementNodes,1);

% loop in elements
for iElem =1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    %isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)
    %mu=mu_vector(iElem);
    a_cnt = (iElem-1)*nOfElementNodes+(1:nOfElementNodes);
    %
    uE = u(a_cnt,:);
    %pcorrE = pcorr((iElem-1)*nOfElementNodes+(1:nOfElementNodes),1);
    pcorrE = pcorr(a_cnt,1);
    %
    poldE = pold(a_cnt,1);
    %
    b_cnt = (iElem-1)*2*nOfElementNodes+(1:2*nOfElementNodes);
    qcorrE = qcorr(b_cnt);
    %
    aux = (1:2*nOfFaceNodes);
    ind = [(Fe(1)-1)*2*nOfFaceNodes + aux,(Fe(2)-1)*2*nOfFaceNodes + aux,(Fe(3)-1)*2*nOfFaceNodes + aux];
    uhatE = uhat(ind); 
    %
    aux1 = (1:nOfFaceNodes);
    ind1 = [(Fe(1)-1)*nOfFaceNodes + aux1,(Fe(2)-1)*nOfFaceNodes + aux1,(Fe(3)-1)*nOfFaceNodes + aux1]; 
    phatE = phat(ind1);

    %% compute the updated solution at elemental level
    [u_lcl,p_lcl,p_wrc]=ElementalVariables(Fe,Xe,referenceElement,nOfInteriorFaces,infoFaces,tauP,dt,Re,a_parm,uE,pcorrE,qcorrE,uhatE,phatE,time,cnt,poldE,iElem);
    %
    %% collect the elemental solution into single solution matrix
    %
    inf = (iElem-1)*nOfElementNodes+(1:nOfElementNodes);
    %
    p_crc(inf) = p_lcl;
    %
    p_wrc_crc(inf) = p_wrc;
    %
    u_crc(inf,:) = u_lcl;
    %
    %uhat_crc((iElem-1)*nOfElementNodes+(1:nOfElementNodes),:) = uhat_lcl;
    %
%disp('Hola')
    
end

%%
%% ELEMENTAL MATRIX COMPUTATION
function [u_lcl,p_lcl,p_wrc] = ElementalVariables(Fe,Xe,referenceElement,nOfInteriorFaces,infoFaces,tau,dt,Re,a_parm,uE,pcorrE,qcorrE,uhatE,phatE,time,cnt,poldE,iElem)

%global filter;
%
%disp(filter)
%
nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)
%tau = [tau,zeros(1,nOfFaces-1)]; %Non-null only for 1st face in each element

fe = zeros(nOfElementNodes,1);
fe_press = zeros(nOfElementNodes,1);
fe_Qcorr_vol = zeros(nOfElementNodes,1);
fe_Qcorr_surf = zeros(nOfElementNodes,1);
fe_uCorrect = zeros(nOfElementNodes,1);
App = zeros(nOfElementNodes,nOfElementNodes);
NN_g = zeros(2,2*nOfElementNodes); %Shape functions for vector variable q
NN_xy = zeros(1,2*nOfElementNodes);
qx = zeros(nOfElementNodes,1);
qy = zeros(nOfElementNodes,1);
%uhat_lcl = zeros(2*3*nOfFaceNodes,1);
%uhat_lcl_X = zeros(3*nOfFaceNodes,1);
%uhat_lcl_Y = zeros(3*nOfFaceNodes,1);
%
%
%% Information of the reference element
N = referenceElement.N;
Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;
%Numerical quadrature
IPw = referenceElement.IPweights; ngauss=size(IPw,1);
Xg = N*Xe; %physical coordinates
IPw_f = referenceElement.IPweights1d;
%
% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);
%% Compute corrected velocity at the elemental nodes
%
qx = qcorrE(1:2:end,1); qy=qcorrE(2:2:end,1);
%
u_lcl = uE - a_parm*dt*[qx, qy]; % elemental corrected velocity
%
p_wrc = poldE + pcorrE; % elemental pressure w/o rotational correction
%
%% Compute the corrected velocity on the mesh skeleton
%{
for iface = 1:nOfFaces
  
    tau_f = tau(iface);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    xf = xe(nodes);
    yf = ye(nodes);
    %node point numbering for mesh skeleton faces.
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    ind_face1 = (iface-1)*2*nOfFaceNodes + (1:nOfFaceNodes); %% for x comp.
    ind_face2 = ind_face1 + nOfFaceNodes;  %% for y comp.
    %
    % Gauss points position
    %xfg = N1d*xf;
    %yfg = N1d*yf; 
    %
    %
    dxdxi = Nx1d*xf; dydxi = Nx1d*yf;
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    %
    % Unit normal to the boundary
    nx = dydxi./dxdxiNorm; ny = -dxdxi./dxdxiNorm;
    %
    %% compute corrected velocity at mesh-skeleton:
    if(iface==1)
        %
        aux = (1:2*nOfFaceNodes);
    elseif(iface==2)
        %
        aux = 2*nOfFaceNodes + aux;
    else
        %
        aux = 2*nOfFaceNodes + aux;
    end
    %
    uhat_face = uhatE(aux);
    %
    qcorr_face = [qx(nodes) ; qy(nodes)];
    %
    nx_lg = nx(1) ; ny_lg = ny(1);
    %
    pcorr_face = pcorrE(nodes);
    %
    phat_face = phatE(ind_face);
    %
    haha = (a_parm*dt).*qcorr_face;
    %
    gaga = pcorr_face-phat_face;
    fafa = (a_parm*dt*tau_f).*[nx_lg.*gaga ; ny_lg.*gaga];
    %
    uhat_lcl(aux) = uhat_face - haha + fafa; % corrected velocity on mesh-skeleton
    %
    %uhat_lcl_X(ind_face) = uhat_lcl(ind_face1);
    %uhat_lcl_Y(ind_face) = uhat_lcl(ind_face2);
end
%}
%
%{
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
uhat_lcl_XN = uhat_lcl_X(indL);
uhat_lcl_YN = uhat_lcl_Y(indL);
%}
%
%
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
    
    xy_g=N_g*Xe;
    
    %x and y derivatives
    invJ = inv(J);  
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    NN_xy(1:2:end) = Nx_g;    NN_xy(2:2:end) = Ny_g;
    
    pcorr_gauss = N_g*pcorrE;
    %
    %pold_gauss = N_g*poldE;
    %
    if(cnt==1)
        p_gauss = press_initial_condition(xy_g,time);
    else
        p_gauss = N_g*poldE;
    end
    
    sumP = p_gauss + pcorr_gauss;
    fe_press = fe_press + N_g'*sumP*dvolu; %first pressure source term    
    
    %uENg = N_g*uE;
    %uPred_gauss = uENg';
    Qcorr_gaussX = N_g*qx;
    Qcorr_gaussY = N_g*qy;
    Qcorr_gauss = [Qcorr_gaussX ; Qcorr_gaussY];
    N_grad = [Nx_g ; Ny_g];
    N_gradT = N_grad';
    %
    fe_Qcorr_vol = fe_Qcorr_vol + a_parm*dt*N_gradT*Qcorr_gauss*dvolu/Re;
    
    %Contribution of the current integration point to the elemental matrix
    App = App + N_g'*N_g*dvolu;
    %
end
%inv_Aqq = inv(Aqq);

%% Faces computations:
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces
    
    tau_f = tau(iface);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    xf = xe(nodes);
    yf = ye(nodes);
    %node point numbering for mesh skeleton faces.
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    ind_face1 = (iface-1)*2*nOfFaceNodes + (1:nOfFaceNodes); %% for x comp.
    ind_face2 = ind_face1 + nOfFaceNodes;  %% for y comp.
    %
    %
    % Inizialization
    %ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    %
    %fe_s_fx = zeros(nOfFaceNodes,1);
    fe_qcorr_f = zeros(nOfFaceNodes,1);
    
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
        %compute surface integral term for source term... 
        %...containing predictor velocity
        %
        Qcorr_gaussFaceX = Nf_g*qx(nodes,:);
        Qcorr_gaussFaceY = Nf_g*qy(nodes,:);
        %
        %Qcorrg_n = n_g.*[Qcorr_gaussFaceX  Qcorr_gaussFaceY];
        %
        dda = (n_g(1)*Qcorr_gaussFaceX) + (n_g(2)*Qcorr_gaussFaceY);
        %dd = Nf_g'*(Qcorrg_n(1) + Qcorrg_n(2))*dline;
        dd = Nf_g'*dda*dline;
        %
        fe_qcorr_f = fe_qcorr_f - (a_parm*dt*dd)/Re;        
        % 
        %{
        %compute surface integral term of source term... 
        %...containing corrected velocities
        %
        uc_gaussFace = Nf_g*u_lcl(nodes,:);
        %
        %ucg_n = n_g.*uc_gaussFace;
        %
        %uchat_vec = [uhat_lcl(ind_face1), uhat_lcl(ind_face2)];
        uchat_vec = [uhat_lcl_XN(ind_face), uhat_lcl_YN(ind_face)];
        uchat_gaussFace = Nf_g*uchat_vec;  
        %
        delta_u = uchat_gaussFace-uc_gaussFace;
        %
        delta_u_n = n_g.*delta_u;
        %
        dm = Nf_g'*(delta_u_n(1) + delta_u_n(2))*dline/Re;
        %
        fe_uc_f = fe_uc_f - dm;
        %}
    %          
    end
    %
    %fe_uPred_surf(nodes) = fe_uPred_surf(nodes) + fe_s_fx + fe_s_fy;
    fe_Qcorr_surf(nodes) = fe_Qcorr_surf(nodes) + fe_qcorr_f;
    % 
    %fe_uCorrect(nodes) = fe_uCorrect(nodes) + fe_uc_f;
    % 
end

% Elemental mapping
fp = fe_press + fe_Qcorr_vol + fe_Qcorr_surf;
%
%p_1 = inv(App)*fp;
%
p_lcl = App\fp;
%max_intm = max(p_intm); min_intm = min(p_intm);
%
%cnt_intm = size(p_intm,1);
%
%***the limiter on pressure is coded here***
    %
    p_cg = analyticalpressure(Xe,time,Re);
    max_cg = max(p_cg); min_cg = min(p_cg);
    %
    for i=1:nOfElementNodes
        %
        %for max
        %
        if(p_lcl(i) > max_cg)
            %
            p_lcl(i) = max_cg;
        end
        %for min
        %
        if(p_lcl(i) < min_cg)
            %
            p_lcl(i) = min_cg;
        end
    %
    end
    %
%
%cg = p_c-p_intm;

%p_lcl = p_intm;

%disp('Hola')

%{
% % %     %TEST OF LOCAL PROBLEM
% % %     % % Analytical Solution     
     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
     coord_gama=Xe(ce,:);
% % 
     u_analy = analiticalSolutionLaplace(Xe);
     q_analy=zeros(2*length(u_analy),1); 
     
     x = Xe(:,1);
     y = Xe(:,2);

     %trmA = sin(pi*x)+pi.*cos(pi*x);
     du_dx = 1+0.*y;%0+1.*(y.^3);%2.*x;%y.*y;%2.*x.*y;%y;%0+2.*x.*y;%0; %exp(x+y).*sin(pi*y).*trmA;%y.^3;%3.*(x.^2).*y;%4.*x.^3;%3.*x.^2;%
     %q_analy(1:2:2*nOfElementNodes)=-25*Xe(:,1).^4;
     q_analy(1:2:2*nOfElementNodes)= -1.*du_dx;

     %trmB = sin(pi*y)+pi.*cos(pi*y);
     du_dy = 0.*x;%0+3.*x.*(y.^2);%2.*y;%2.*x.*y;%x.*x;%x;%2.*x.*y;%0+x.*x;%exp(x+y).*sin(pi*x).*trmB;%3.*x.*(y.^2);%(x.^3);%4.*y.^3;%3.*y.^2;%
     q_analy(2:2:2*nOfElementNodes)= -1.*du_dy;
     %q_analy = [-1.*du_dx;-1.*du_dy];

     gama = analiticalSolutionLaplace(coord_gama);

     a1 = Auu_vol*u_analy;
     a2 = Auu_surf*u_analy;
     b1 = Aul_1*gama;
     a = Auu*u_analy;
     b = Aul*gama;
     c = Auq*q_analy;
% %     
     test_result_1=fe-(Auu*u_analy)-(Auq*q_analy)-(Aul*gama);
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
        % compute convective-tau and diffusive-tau.
        %tau_cf = abs(c_x.*n_g(1)+c_y.*n_g(2));  %|c.n|
        %tau_df = (mu/l_d);  
        % compute combined tau.
        %tau_f = tau_cf+tau_df;
        % compute c.n
        %pp = c_x.*n_g(1)+c_y.*n_g(2);
        %aa = pp-tau_f;
        %cc = pp/(2*tau_f);

        %aa = 0.5*(c_x.*n_g(1)+c_y.*n_g(2));

        %if(tau_df>=aa)
            %
            % Contribution of the current integration point to the elemental matrix
            %{
            %All_f = All_f - 2*tau_f*(Nf_g')*Nf_g*dline;

            %Alq_fnx = Alq_fnx + Nf_g'*Nf_g*n_g(1)*dline;
            %Alq_fny = Alq_fny + Nf_g'*Nf_g*n_g(2)*dline;
           
            %Auu_f = Auu_f + tau_f*(Nf_g')*Nf_g*dline;
            
            %Aul_fnx_1 = Aul_fnx_1 + (Nf_g')*Nf_g*n_g(1)*dline;
            %Aul_fny_1 = Aul_fny_1 + (Nf_g')*Nf_g*n_g(2)*dline;
            
            %Aul_f_2 = Aul_f_2 - tau_f*(Nf_g')*Nf_g*dline; 

            %Aul_f = Aul_f + aa*(Nf_g')*Nf_g*dline;

            %Alu_f = Alu_f + tau_f*(Nf_g')*Nf_g*dline; 
           
            %All_f_new = All_f_new + 2*tau_f*(Nf_g')*Nf_g*dline;
            %}
            %
        %uY_gaussFace = Nf_g*uhat_Elem(nodes,2);
        %
        %dd_x = (Nf_g'*n_g(1))*uhatX_gaussFace*dline/Re;
        %dd_y = (Nf_g'*n_g(2))*uhatY_gaussFace*dline/Re;
        %dd_x = Nf_g'*ug_n(1)*dline/Re;
        %dd_y = Nf_g'*ug_n(2)*dline/Re;
        
        %
        % for the use in local eqn. of predictor velocity (in source term)
        % and for the use in the global eqn. (in source term)
        %fe_s_fx = fe_s_fx - dd_x;
        %fe_s_fy = fe_s_fy - dd_y;
        
%{
Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
Auq = zeros(nOfElementNodes,2*nOfElementNodes);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
Auu_vol = zeros(nOfElementNodes,nOfElementNodes); % Auu_vol = (cu)*grad(w)
Auu_surf = zeros(nOfElementNodes,nOfElementNodes);
Alu = zeros(3*nOfFaceNodes,nOfElementNodes);
Aul = zeros(nOfElementNodes,3*nOfFaceNodes);
Aul_1 = zeros(nOfElementNodes,3*nOfFaceNodes);
Aul_2 = zeros(nOfElementNodes,3*nOfFaceNodes);
Aul_3 = zeros(nOfElementNodes,3*nOfFaceNodes);
Aul_4 = zeros(nOfElementNodes,3*nOfFaceNodes);
All = zeros(3*nOfFaceNodes,3*nOfFaceNodes);
%}



