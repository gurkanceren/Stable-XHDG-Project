function [Auq,Aqu,Alu,Alq,Auu,Aqq,All,fe]=xhdg_matrix(referenceElement,mu,tau,Xe,LSe,iElem,F)

nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);
faceNodes = referenceElement.faceNodes;

tau = repmat(tau,1,nOfFaces); %Uniform tau parameter (all faces)


% Information of the reference element
N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;
N = referenceElement.N;
Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
IPw = referenceElement.IPweights;
IPw_f = referenceElement.IPweights1d;
ngauss=size(IPw,1);
weight=IPw;

NN_g = zeros(2,2*nOfElementNodes); %Shape functions for vector variable q
NN_xy = zeros(1,2*nOfElementNodes);


% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

% coordinates of gauss pts on real element
Xg = N*Xe ;

% analytical laplacian
sourceTerm = @analiticalLaplacianLaplace;


% Initialization

fe  = zeros(nOfElementNodes,1);
Aqq = zeros(2*nOfElementNodes,2*nOfElementNodes);
Auq = zeros(nOfElementNodes,2*nOfElementNodes);
Aqu = zeros(2*nOfElementNodes,nOfElementNodes);
Alq = zeros(3*nOfFaceNodes,2*nOfElementNodes);
Auu = zeros(nOfElementNodes,nOfElementNodes);
Alu = zeros(3*nOfFaceNodes,nOfElementNodes);
All = zeros(3*nOfFaceNodes,3*nOfFaceNodes);

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
    
    %Contribution of the current integration point to the elemental matrix
    Aqq = Aqq + NN_g'*NN_g*dvolu/mu;
    Auq = Auq + N_g' *NN_xy*dvolu;
    Aqu = Aqu - NN_xy'*N_g*dvolu;
    fe = fe + N_g'*sourceTerm(xy_g)*-mu*dvolu; 
end


%% Faces computations:
ngauss_f = length(IPw_f);
for iface = 1:nOfFaces
        
    tau_f = tau(iface);
    
    % Nodes in the face
    nodes = faceNodes(iface,:);
    xf = xe(nodes);
    yf = ye(nodes);
    
    % Inizialization
    ind_face = (iface-1)*nOfFaceNodes + (1:nOfFaceNodes);
    Alq_fnx = zeros(nOfFaceNodes,nOfFaceNodes);
    Alq_fny = zeros(nOfFaceNodes,nOfFaceNodes);
    Auu_f = zeros(nOfFaceNodes,nOfFaceNodes);
    All_f = zeros(nOfFaceNodes);
    
        
        glim=ngauss_f;
        %Shape functions and derivatives
        Nf_g = N1d;
        Nfxi_g = Nx1d;
        %Integration weight
        IPwFace=IPw_f; 
       
    
    for g = 1:glim
        
        % Integration weight
        xyDer_g = Nfxi_g(g,:)*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline=IPwFace(g)*xyDerNorm_g;
        
        % Unit normal to the boundary
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) -t_g(1)];
        
        % Contribution of the current integration point to the elemental matrix
        Alq_fnx = Alq_fnx + Nf_g(g,:)'*Nf_g(g,:)*n_g(1)*dline;
        Alq_fny = Alq_fny + Nf_g(g,:)'*Nf_g(g,:)*n_g(2)*dline;
        Auu_f = Auu_f + tau_f *mu*Nf_g(g,:)'*Nf_g(g,:)*dline;
        All_f = All_f - 2*mu*tau_f *Nf_g(g,:)'*Nf_g(g,:)*dline;
    end
    
    
    Alq(ind_face,2*nodes-1) = Alq(ind_face,2*nodes-1) + Alq_fnx;
    Alq(ind_face,2*nodes) = Alq(ind_face,2*nodes) + Alq_fny;
    Auu(nodes,nodes) = Auu(nodes,nodes) + Auu_f;
    Alu(ind_face,nodes) = Alu(ind_face,nodes) + Auu_f;
    All(ind_face,ind_face) = All(ind_face,ind_face) + All_f;
    
end
     Aul = -Alu'; Aql = Alq';
    
% %     
% % %     %TEST OF LOCAL PROBLEM
% % %     % % Analytical Solution     
% %     ce=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
% %     coord_gama=Xe(ce,:);
% % 
% %     u_analy = analiticalSolutionLaplace(Xe);
% %     q_analy=zeros(2*length(u_analy),1); q_analy(1:2:2*nOfElementNodes)=-240*Xe(:,1).^5;
% %     gama = analiticalSolutionLaplace(coord_gama);
% %     
% %     test_result_1=fe-(Auu*u_analy)-(Auq*q_analy)-(Aul*gama);
% %     test_result_2=(Aqq*q_analy)+(Aqu*u_analy)+(Aql*gama);
% %     
% %     if any([max(abs(test_result_1)),max(abs(test_result_2))]>1.e-14)
% %             disp(sprintf('      LP max(abs(residual)) eq. #1: %e',max(abs(test_result_1))));
% %             disp(sprintf('      LP max(abs(residual)) eq. #2: %e',max(abs(test_result_2))));
% %     else
% %         disp('LP OK')
% %     end
% %     
    
 