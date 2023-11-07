function [Aqqnew,Auqnew,Aqunew,Alu,Alq,fe,All,Auunew]=xhdg_matrix_cut(referenceElement,...
    mu1,mu2,tau,Xe,LSe,iElem,F)


nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;
IPw_f = referenceElement.IPweights1d;
tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)

gama = zeros(6*nOfFaceNodes,1); %for TEST

[Nold,Nxi_old,Neta_old,Nenr,Nxienr,Netaenr,weight,ngauss,mu,PtsInt,FaceInfo]=modifiedReferenceElement(referenceElement,LSe,mu1,mu2);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% physical coordinates
Xg = Nold*Xe ;

% analytical laplacian
sourceTerm = @analiticalLaplacianLaplace;

NN_g = zeros(2,4*nOfElementNodes); %Shape functions for vector variable q
NN_xy = zeros(1,4*nOfElementNodes);

%Initialization
Alq = zeros(6*nOfFaceNodes,4*nOfElementNodes);
Auu = zeros(2*nOfElementNodes,2*nOfElementNodes);
Alu = zeros(6*nOfFaceNodes,2*nOfElementNodes);
All = zeros(6*nOfFaceNodes,6*nOfFaceNodes);
Aqq = zeros(4*nOfElementNodes,4*nOfElementNodes);
Auq = zeros(2*nOfElementNodes,4*nOfElementNodes);
Aqu = zeros(4*nOfElementNodes,2*nOfElementNodes);
fe = zeros(2*nOfElementNodes,1);


%% Volume computation: loop in Gauss points
for g = 1:ngauss
    
    %Shape functions and derivatives at the current integration point
    N_g = Nenr(g,:);  N_gold = Nold(g,:);
    NN_g(1,1:2:end)=N_g; NN_g(2,2:2:end)=N_g;
    Nxi_g = Nxienr(g,:);    Neta_g = Netaenr(g,:);
    Nxi_g_old=Nxi_old(g,:);  Neta_g_old=Neta_old(g,:);
    
    %Jacobian
    J = [Nxi_g_old*xe  Nxi_g_old*ye
        Neta_g_old*xe  Neta_g_old*ye];
    
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end
    
    %Integration weight
    dvolu=weight(g)*det(J);
    
    %integration pt coord
    xy_g=N_gold*Xe;
    
    %x and y derivatives
    invJ = inv(J);  %Attention!!!!!
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    NN_xy(1:2:end) = Nx_g;    NN_xy(2:2:end) = Ny_g;
    
    %Contribution of the current integration point to the elemental matrix
    
    Aqq = Aqq + NN_g'*NN_g*dvolu/mu(g);
    Auq = Auq + N_g'*NN_xy*dvolu;
    Aqu = Aqu - NN_xy'*N_g*dvolu;
    fe = fe + N_g'*sourceTerm(xy_g)*-mu(g)*dvolu;      % !!!!source term here !!!!!
end



%% Faces computations:
ngauss_f = length(IPw_f);
N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;

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
    ind_face =(iface-1)*(2*nOfFaceNodes) + (1:(2*nOfFaceNodes));
    Alq_fnx = zeros(2*nOfFaceNodes,2*nOfFaceNodes);
    Alq_fny = zeros(2*nOfFaceNodes,2*nOfFaceNodes);
    Auu_f = zeros(2*nOfFaceNodes,2*nOfFaceNodes);
    Alu_f = zeros(2*nOfFaceNodes,2*nOfFaceNodes);
    All_f = zeros(2*nOfFaceNodes,2*nOfFaceNodes);
    
    if isempty(k) && any(LSe(faceNodes(iface,:))>0) ;  % NORMAL FACES in D1
        
        glim=ngauss_f;
        
        %Shape functions and derivatives
        Nf_g = [N1d N1d*0];
        Nenr_g = [N1d N1d*-1];
        Nfxi_g_old = Nx1d;
        %Integration weight
        IPwFace=IPw_f;
        mu=ones((ngauss_f),1)*mu1;
        %gama(ind_face) = [analiticalSolutionLaplace([xf,yf]);zeros(nOfFaceNodes,1)]; %for TEST
        
    elseif  isempty(k) && any(LSe(faceNodes(iface,:))<0) ; %Normal face in D2
        
        glim=ngauss_f;
        
        %Shape functions and derivatives
        Nf_g = [N1d N1d*0];
        Nenr_g = [N1d N1d*1];
        Nfxi_g_old = Nx1d;
        %Integration weight
        IPwFace=IPw_f;
        mu=ones((ngauss_f),1)*mu2;
        %gama(ind_face) = [analiticalSolutionLaplace([xf,yf]);zeros(nOfFaceNodes,1)]; %for TEST
        
    else   ~isempty(k);    %CUT FACE
        
        % Enriched Shape Functions for cut faces:
        
        [Nx_old,Nx_enrf,N_enrf,wgp_f,glim]=modifiedReferenceElement1Cut(referenceElement,LSe,nodes);
        
        Nf_g= N_enrf;
        Nenr_g = N_enrf;
        Nfxi_g_old=Nx_old;
        %Integration weight
        IPwFace=wgp_f;
        mul=ones((ngauss_f),1)*mu1;
        mur=ones((ngauss_f),1)*mu2;
        mu=[mul;mur];
        %gama(ind_face) = analyticalSolutionWithHeaviside([xf,yf]); %for TEST
    end
    
    %  LOOP IN GAUSS POINTS
    for g = 1:glim
        
        % Integration weight
        xyDer_g = Nfxi_g_old(g,:)*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline=IPwFace(g)*xyDerNorm_g;
        
        % Unit normal to the boundary
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) -t_g(1)];
        
        % Contribution of the current integration point to the elemental matrix
        Alq_fnx = Alq_fnx + Nf_g(g,:)'*Nenr_g(g,:)*n_g(1)*dline;
        Alq_fny = Alq_fny + Nf_g(g,:)'*Nenr_g(g,:)*n_g(2)*dline;
        Auu_f = Auu_f + tau_f*mu(g)*Nenr_g(g,:)'*Nenr_g(g,:)*dline;
        Alu_f = Alu_f + tau_f*mu(g)*Nf_g(g,:)'*Nenr_g(g,:)*dline;
        All_f = All_f - tau_f*2*mu(g)*Nf_g(g,:)'*Nf_g(g,:)*dline;
        
    end
    
    nodes = [nodes nodes+nOfElementNodes];
    
    Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Alq_fnx;
    Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Alq_fny;
    Auu(nodes,nodes) = Auu(nodes,nodes) + Auu_f;
    Alu(ind_face,nodes) = Alu(ind_face,nodes) + Alu_f;
    All(ind_face,ind_face) = All(ind_face,ind_face) + All_f;
    
end
    
    %% INTEGRATIONS OVER I
    
    
    %Calculation of new Aquutilda,Auutilda and AuuonI matrices because of applied
    %Neumann boundary conditions on I
    
    p=referenceElement.degree;
    g=length(referenceElement.IPweights1d);    
    PtsInt=reshape(PtsInt',2,p+1)'; %interface p+1 nodes in the REFERENCE element
    zg=referenceElement.N1d*PtsInt;%integration points on the interface in the REFERENCE element
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zg);
    Nzg = shapeFunctions(:,:,1)';  % 2d shape functions on the refence element at integration points
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,PtsInt);
    NPtsInt = shapeFunctions(:,:,1)';  %2D Shape functions on the REFERENCE element at interface nodes
    
    %Enriched Shape Functions:
    
    NzgR=[Nzg Nzg*1];
    NzgL=[Nzg Nzg*-1];
    
    Nref=referenceElement.N1d;
    
    muL=ones(g,1)*mu1;
    muR=ones(g,1)*mu2;
    
    Pphy=NPtsInt*Xe; %p+1 interface nodes on the Physical element
    Iprime=referenceElement.N1dxi*Pphy;
    
    %Initialization
    
    AuuonI=zeros(2*nOfElementNodes,2*nOfElementNodes);
    Auutilda=zeros(2*nOfElementNodes,nOfFaceNodes);
    Aqutilda=zeros(4*nOfElementNodes,nOfFaceNodes);
    Autildaq=zeros(nOfFaceNodes,4*nOfElementNodes);
    Aqutilda_x_f=zeros(2*nOfElementNodes,nOfFaceNodes);
    Aqutilda_y_f=zeros(2*nOfElementNodes,nOfFaceNodes);
    Autildaq_x_f=zeros(nOfFaceNodes,2*nOfElementNodes);
    Autildaq_y_f=zeros(nOfFaceNodes,2*nOfElementNodes);
    Autildautilda=zeros(nOfFaceNodes,nOfFaceNodes);
    
    for igauss=1:g
        
        
        % Integration weight
        normIprime = norm(Iprime(igauss,:));
        wg=referenceElement.IPweights1d(igauss)*normIprime;
        
        % Unit normal to the boundary
        t_g = Iprime(igauss,:)/normIprime;
        n_g = [t_g(2) -t_g(1)];
        
        
        AuuonI= AuuonI+tau_f*(NzgR(igauss,:)'*NzgR(igauss,:)*muR(igauss)+NzgL(igauss,:)'*NzgL(igauss,:)*muL(igauss))*wg;
        Auutilda= Auutilda-tau_f*(NzgR(igauss,:)'*Nref(igauss,:)*muR(igauss)+NzgL(igauss,:)'*Nref(igauss,:)*muL(igauss))*wg;
        Aqutilda_x_f =  Aqutilda_x_f+( (NzgL(igauss,:)-NzgR(igauss,:))'*Nref(igauss,:) )*(n_g(1)*wg);
        Aqutilda_y_f =  Aqutilda_y_f+( (NzgL(igauss,:)-NzgR(igauss,:))'*Nref(igauss,:) )*(n_g(2)*wg);
        Autildautilda= Autildautilda-tau_f*Nref(igauss,:)'*Nref(igauss,:)*(muR(igauss)+muL(igauss))*wg;
         
         Autildaq_x_f =  Autildaq_x_f+(Nref(igauss,:)'*NzgL(igauss,:)*n_g(1) - Nref(igauss,:)'*NzgR(igauss,:)*n_g(1))*wg;
         Autildaq_y_f =  Autildaq_y_f+(Nref(igauss,:)'*NzgL(igauss,:)*n_g(2) - Nref(igauss,:)'*NzgR(igauss,:)*n_g(2))*wg;
        
    end
    
    Aqutilda(1:2:end,:) = Aqutilda_x_f;
    Aqutilda(2:2:end,:) = Aqutilda_y_f;
    
    Autildaq(:,1:2:end) = Autildaq_x_f;
    Autildaq(:,2:2:end) = Autildaq_y_f;
    Autildau=-Auutilda';
    
    Me=Autildautilda^-1;
    Bq=-Me*Autildaq;
    Bu=-Me*Autildau;
    
    
    Auqnew=Auq+Auutilda*Bq;
    Auunew=Auu+AuuonI+Auutilda*Bu;
    Aqunew=Aqu+Aqutilda*Bu;
    Aqqnew=Aqq+Aqutilda*Bq;
    
    
    
     % %     %smalltest1
% %     
% %     
% %     utilda=(40*0.1875^6)*ones(nOfFaceNodes,1);
% %     u_analy=analyticalSolutionWithHeaviside(Xe);
% %     
% %     q_analy=zeros(2*length(u_analy),1);
% %     q_analy(1:2:2*nOfElementNodes)=-240*Xe(:,1).^5; 
% %     
% %     % % Test for two equations in local problem 
% %     
% %     new_test_eqn1=[Auu+AuuonI]*u_analy + Auq*q_analy+Aul*gama+Auutilda*utilda-fe;  % satisfied
% %     new_test_eqn2=Aqu*u_analy + Aqq*q_analy+Aql*gama+Aqutilda*utilda;              % satisfied
% %     
% %     % % Test for conservativity condition on I 
% %     result=Autildau*u_analy+Autildaq*q_analy+Autildautilda*utilda;   %%% satisfied  !!
% %    
% %      % % 
% %     
% %     test_result_1=fe-(Auunew*u_analy)-(Auqnew*q_analy)-(Aul*gama);  %% satisfied
% %     test_result_2= (Aqqnew*q_analy)+(Aqunew*u_analy)+(Aql*gama);    %% satisfied 
% % 
% %     if any([max(abs(test_result_1)),max(abs(test_result_2))]>1.e-14)
% %             disp(sprintf('(cut) LP max(abs(residual)) eq. #1: %e',max(abs(test_result_1))));
% %             disp(sprintf('      LP max(abs(residual)) eq. #2: %e',max(abs(test_result_2))));
% %     else
% %         disp('LP OK')
% %     end
    
  