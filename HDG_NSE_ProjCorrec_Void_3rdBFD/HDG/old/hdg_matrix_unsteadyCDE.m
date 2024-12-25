function [KK,f, QQ, UU, Qf, Uf, uDirichlet] = hdg_matrix(X,T,F,referenceElement,infoFaces,tau,mu_vector,upwind,dt,cnt,u,time,c_x,c_y)
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

%u_old = ProjectingOldSolution(u,nOfElements,X,T,referenceElement);

% loop in elements
for iElem =1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)
    mu=mu_vector(iElem);
    u_Elem = u((iElem-1)*nOfElementNodes+(1:nOfElementNodes));

    % elemental matrices
    [Qe,Ue,Qfe,Ufe,Alq,Alu,All] = KKeElementalMatrices(Xe,referenceElement,tau,mu,upwind,dt,cnt,u_Elem,time,c_x,c_y);
    
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
    KKe = Alq*Qe + Alu*Ue + 0.5*All;
    ffe = -(Alq*Qfe + Alu*Ufe);
    
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
uDirichlet=computeProjectionFacesDirchlet(infoFaces.extFaces,X,T,referenceElement,time);
% System reduction (Dirichlet faces are set to prescribed value)
nDOF = nOfInteriorFaces*nOfFaceNodes;
f = f(1:nDOF)-KK(1:nDOF,nDOF+1:end)*uDirichlet;
KK=KK(1:nDOF,1:nDOF);

%%
%% ELEMENTAL MATRIX
function [Q,U,Qf,Uf,Alq,Alu,All] = KKeElementalMatrices(Xe,referenceElement,tau,mu,upwind,dt,cnt,u_Elem,time,c_x,c_y)
                                                        
nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)
%tau = [tau,zeros(1,nOfFaces-1)]; %Non-null only for 1st face in each element

fe_s = zeros(nOfElementNodes,1);
fe_uns = zeros(nOfElementNodes,1);
u_old = zeros(3,1);
Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
Auq = zeros(nOfElementNodes,2*nOfElementNodes);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
Auu_vol_s = zeros(nOfElementNodes,nOfElementNodes); % Auu_vol = (cu)*grad(w)
Auu_vol_t = zeros(nOfElementNodes,nOfElementNodes);
Auu_surf = zeros(nOfElementNodes,nOfElementNodes);
Alu = zeros(3*nOfFaceNodes,nOfElementNodes);
Aul = zeros(nOfElementNodes,3*nOfFaceNodes);
Aul_1 = zeros(nOfElementNodes,3*nOfFaceNodes);
Aul_2 = zeros(nOfElementNodes,3*nOfFaceNodes);
Aul_3 = zeros(nOfElementNodes,3*nOfFaceNodes);
Aul_4 = zeros(nOfElementNodes,3*nOfFaceNodes);
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

l_d = 1; % representative diffusive length scale
%c_x = 1; % advection speed of the signal in x-direction
%c_y = 1; % advection speed of the signal in y-direction

% call the subroutine which computes the source term.
sourceTerm = @analiticalLaplacianLaplace;

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

    % advection speed:
    %c_x = -4.*xy_g(:,2);
    %c_y = 4.*xy_g(:,1);
    
    %x and y derivatives
    invJ = inv(J);  
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    NN_xy(1:2:end) = Nx_g;    NN_xy(2:2:end) = Ny_g;
    
    %Contribution of the current integration point to the elemental matrix
    Aqq = Aqq + NN_g'*NN_g*dvolu/mu;
    Auq = Auq + N_g'*NN_xy*dvolu;
    %Auu_vol_x = N_g'*Nx_g*dvolu;
    %Auu_vol_y = N_g'*Ny_g*dvolu;

    %Auu_vol_x1 = N_g'*c_x*Nx_g*dvolu;
    %Auu_vol_y1 = N_g'*c_y*Ny_g*dvolu;
    %Auu_vol_x = -c_x*N_g'*Nx_g*dvolu;
    %Auu_vol_y = -c_y*N_g'*Ny_g*dvolu;
    c_vec = [c_x c_y];
    cu = c_vec'*N_g;
    %cu_t = cu';

    grad_w = [Nx_g ; Ny_g];
    grad_w_t = grad_w';

    temp = -grad_w_t*cu;

    Auu_vol_s = Auu_vol_s + temp*dvolu;
    Auu_vol_t = Auu_vol_t + (N_g'*N_g*dvolu/dt);
    % fe = fe - N_g'*sourceTerm(xy_g)*mu*dvolu;
    %ha = sourceTerm(xy_g,c_x,c_y,time);
    fe_s = fe_s + N_g'*sourceTerm(xy_g,c_x,c_y,time)*dvolu;
    
    if(cnt==1)
        u_old_gauss = initial_condition(xy_g);
        %u_old_gauss = analiticalSolutionLaplace(xy_g,time);
    else
        u_old_gauss = N_g*u_Elem; 
    end 
    %u_old_gauss = N_g*u_Elem;
   
    fe_uns = fe_uns + N_g'*u_old_gauss*dvolu/dt;
end
%inv_Aqq = inv(Aqq);

%% Faces computations:
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces
    
    % tau_f = tau(iface);
    
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
    Alu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_f = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_fnx_1 = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_fny_1 = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_f_2 = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_fnxx = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_fnxy = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_fnyx = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_fnyy = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_fnx_4 = zeros(nOfFaceNodes,nOfFaceNodes);
    Aul_fny_4 = zeros(nOfFaceNodes,nOfFaceNodes);
    Auu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    %All_f = zeros(nOfFaceNodes);
    All_f = zeros(nOfFaceNodes,nOfFaceNodes);
    All_f_new = zeros(nOfFaceNodes,nOfFaceNodes);
    
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

        % compute convective and diffusive eta:
        eta_c = abs(c_x.*n_g(1)+c_y.*n_g(2));  % |c.n|
        eta_d = (mu/l_d);
        
        
        if (upwind==1)
            % compute c.n+ and |c.n+|:
            %disp('Hola')
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
            %disp('Lola')
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
            % Auu_f = Auu_f + tau_f *mu*(Nf_g')*Nf_g*dline;
            % All_f = All_f - tau_f *2*mu*(Nf_g')*Nf_g*dline;
            All_f = All_f - tau_f *2*(Nf_g')*Nf_g*dline;

            Alq_fnx = Alq_fnx + Nf_g'*Nf_g*n_g(1)*dline;
            Alq_fny = Alq_fny + Nf_g'*Nf_g*n_g(2)*dline;

            %Alu_f = Alu_f + tau_f*(Nf_g')*Nf_g*dline; 

            Auu_f = Auu_f + tau_f*(Nf_g')*Nf_g*dline;

            Aul_f = Aul_f + bb*(Nf_g)'*Nf_g*dline;
            
            %Aul_fnx_1 = Aul_fnx_1 + (Nf_g')*Nf_g*c_x*n_g(1)*dline;
            %Aul_fny_1 = Aul_fny_1 + (Nf_g')*Nf_g*c_y*n_g(2)*dline;

            %Aul_f_2 = Aul_f_2 - tau_f*(Nf_g')*Nf_g*dline;

            %Aul_fnxx = Aul_fnxx + (0.5/tau_f)*c_x*(Nf_g')*Nf_g*n_g(1)*n_g(1)*dline;
            %Aul_fnxy = Aul_fnxy + (0.5/tau_f)*c_x*(Nf_g')*Nf_g*n_g(1)*n_g(2)*dline;
            %Aul_fnyx = Aul_fnyx + (0.5/tau_f)*c_y*(Nf_g')*Nf_g*n_g(2)*n_g(1)*dline;
            %Aul_fnyy = Aul_fnyy + (0.5/tau_f)*c_y*(Nf_g')*Nf_g*n_g(2)*n_g(2)*dline;

            %Aul_fnx_4 = Aul_fnx_4 - 0.5*(Nf_g')*Nf_g*n_g(1)*dline;
            %Aul_fny_4 = Aul_fny_4 - 0.5*(Nf_g')*Nf_g*n_g(2)*dline;
            
            %Aul_f1 = Aul_f1 + pp*(Nf_g')*Nf_g*dline;

            All_f_new = All_f_new + 2*pp*(Nf_g')*Nf_g*dline;
        else
            disp('Lemma 3.1 violated');
        end
       
    end
    Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Alq_fnx;
    Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Alq_fny;

    %Alu(ind_face,nodes) = Alu(ind_face,nodes) + Alu_f;

    Aul(nodes,ind_face) = Aul(nodes,ind_face) + Aul_f;
    %Aul_1(nodes,ind_face) = Aul_1(nodes,ind_face) + Aul_fnx_1 + Aul_fny_1;
    %Aul_2(nodes,ind_face) = Aul_2(nodes,ind_face) + Aul_f_2;
    %Aul_3(nodes,ind_face) = Aul_3(nodes,ind_face)+Aul_fnxx+Aul_fnxy+Aul_fnyx+Aul_fnyy;
    %Aul_4(nodes,ind_face) = Aul_4(nodes,ind_face)+Aul_fnx_4+Aul_fny_4;
    
    Auu_surf(nodes,nodes) = Auu_surf(nodes,nodes) + Auu_f;
    Alu(ind_face,nodes) = Alu(ind_face,nodes) + Auu_f;
    All(ind_face,ind_face) = All(ind_face,ind_face) + All_f + All_f_new;
end


% Elemental mapping
Aqu = -Auq'; Aql = Alq'; %Aul = -Alu'; 
%Aul = Aul_1 + Aul_2;
Auu = Auu_vol_t + Auu_vol_s + Auu_surf;


A = [Auu Auq; Aqu Aqq];
UQ = -A\[Aul;Aql];

fq = zeros(2*nOfElementNodes,1);
fUQ= A\[(fe_uns+fe_s);fq];
%fUQ= A\[(fe_s);fq];

U = UQ(1:nOfElementNodes,:); Uf=fUQ(1:nOfElementNodes); % maps lamba into U
Q = UQ(nOfElementNodes+1:end,:); Qf=fUQ(nOfElementNodes+1:end); % maps lamba into Q 



%
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






