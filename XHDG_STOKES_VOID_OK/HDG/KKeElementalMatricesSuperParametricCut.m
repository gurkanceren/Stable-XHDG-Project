function [LocalSolverMatElem,LocalSolverVecElem,Alu,AlL,Alp,All,Arl] = KKeElementalMatricesSuperParametricCut(mu,Xe,referenceElement,tau,LSe)

nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
faceNodesApprox = referenceElement.faceNodes;
faceNodes = referenceElement.faceNodes;
nOfFaces = 3; %triangles
gama=[];
%Information of the reference element
% N = referenceElement.N;
% Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
% NxiGeo = referenceElement.NxiGeo; NetaGeo = referenceElement.NetaGeo;
% NGeo = referenceElement.NGeo;
% Nx1dGeo = referenceElement.N1dxiGeo;
% IPw = referenceElement.IPweights; ngauss = length(IPw);

[N,Nxi,Neta,weight,ngauss,mu,PtsInt,FaceInfo]=modifiedReferenceElement(referenceElement,LSe,mu);
referenceElement.IPweights=weight';

%%Volume computations
% Jacobian
% J11 = referenceElement.Nxi*Xe(:,1); J12 = referenceElement.Nxi*Xe(:,2); 
% J21 = referenceElement.Neta*Xe(:,1); J22 = referenceElement.Neta*Xe(:,2); 
% detJ = J11.*J22-J12.*J21;

J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); 
J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); 
detJ = J11.*J22-J12.*J21;

%maybe we should use bsxfun instead of diagonal matrices...
%dvolu = spdiags(referenceElement.IPweights.*detJ(1:ngauss),0,ngauss,ngauss);  % !!!
dvolu = spdiags(weight'.*detJ,0,ngauss,ngauss);  % !!!

invJ11 = spdiags(J22./detJ,0,ngauss,ngauss);
invJ12 = spdiags(-J12./detJ,0,ngauss,ngauss);
invJ21 = spdiags(-J21./detJ,0,ngauss,ngauss);
invJ22 = spdiags(J11./detJ,0,ngauss,ngauss);
% xy-derivatives for approximation
Nx = invJ11*Nxi + invJ12*Neta;
Ny = invJ21*Nxi + invJ22*Neta;

%Computation of r.h.s. source term (analytical laplacian)
Xg = N*Xe;
% figure(1); hold on;
% plot(Xg(:,1),Xg(:,2),'*b')
%Xg = referenceElement.N*Xe;
sourceTerm = sourceTermStokes(Xg,mu); %Missing!!!
fe = N'*(dvolu*sourceTerm);
fe = reshape(fe,2*nOfElementNodes,1);
%Basic elemental matrices
Me = N'*(dvolu*N); Cxe=N'*(dvolu*Nx); Cye=N'*(dvolu*Ny); 
%HDG elemental matrices
aux1=1:nOfElementNodes; aux2=aux1+nOfElementNodes; aux3=aux2+nOfElementNodes; aux4=aux3+nOfElementNodes;
ALL = zeros(4*nOfElementNodes,4*nOfElementNodes);
 ALL(aux1,aux1)=Me; ALL(aux2,aux2)=Me; ALL(aux3,aux3)=Me; ALL(aux4,aux4)=Me; 
ALu = zeros(4*nOfElementNodes,2*nOfElementNodes);
 ALu(aux1,aux1)=Cxe'; ALu(aux2,aux1)=Cye'; ALu(aux3,aux2)=Cxe'; ALu(aux4,aux2)=Cye';
AuL = -mu*ALu';
Aup = [Cxe;Cye]; Apu = Aup';


%1D Numerical quadrature
N1d = referenceElement.N1d; Nx1d = referenceElement.Nxi1d;
IPw_f = referenceElement.IPweights1d; ngf = length(IPw_f);

%% Faces computations
ALl = zeros(4*nOfElementNodes,3*2*nOfFaceNodes);
Auu = zeros(2*nOfElementNodes,2*nOfElementNodes);
Aul = zeros(2*nOfElementNodes,3*2*nOfFaceNodes);
Apl = zeros(nOfElementNodes,3*2*nOfFaceNodes);
All = zeros(3*2*nOfFaceNodes,3*2*nOfFaceNodes);
Arp = zeros(1,nOfElementNodes);
Arl = zeros(1,3*2*nOfFaceNodes);


 for iface = 1:nOfFaces
   
    if ~isempty(FaceInfo)
        k=find(FaceInfo(:,1)==iface);
    else
        k=[];
    end
    
    Xf = Xe(faceNodes(iface,:),:); % Nodes in the face
    
    if ~isempty(k);  %Cut face 
        
        %Shape Functions for cut faces:
        [zgp_f,wgp_f,n1_f,n2_f,IntPt] = ModifyQuadrature1D(LSe(referenceElement.faceNodes(iface,:),1),referenceElement);
        shapeFunctions_f=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord1d,zgp_f);        
        ngf=n1_f;
        % Shape functions and derivatives
        N1d =shapeFunctions_f(:,:,1)';    Nx1d = shapeFunctions_f(:,:,2)';
        N1d = N1d(1:n1_f,:);     Nx1d = Nx1d(1:n1_f,:);
        %Integration weight
        IPw_f=wgp_f(1:n1_f);
        dxdxi = Nx1d*Xf(:,1); dydxi = Nx1d*Xf(:,2);
        gama=[gama; analyticalVelocityStokes(Xf,mu)];
        
    elseif isempty(k) && any(LSe(referenceElement.faceNodes(iface,:))>0) ;  %Material face 

           N1d = referenceElement.N1d; Nx1d = referenceElement.Nxi1d;
           IPw_f = referenceElement.IPweights1d; ngf = length(IPw_f);
           dxdxi = Nx1d*Xf(:,1); dydxi = Nx1d*Xf(:,2);
           gama=[gama; analyticalVelocityStokes(Xf,mu)];
        
    elseif  isempty(k) && any(LSe(referenceElement.faceNodes(iface,:))<0) ;  %Void face

        
        IPw_f=IPw_f.*0;
        dxdxi = Nx1d*Xf(:,1); dydxi = Nx1d*Xf(:,2);
        gama=[gama; analyticalVelocityStokes(Xf,mu)];
        
    end
   
    tau_f = tau(iface);
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = dxdxiNorm.*IPw_f;
    nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm;

    %Face basic matrices
    Mf = N1d'*(spdiags(dline,0,ngf,ngf)*N1d)*tau_f;
    Mfnx = N1d'*(spdiags(dline.*nx,0,ngf,ngf)*N1d);
    Mfny = N1d'*(spdiags(dline.*ny,0,ngf,ngf)*N1d);
    intN = dline'*N1d;
    intNnx = (dline.*nx)'*N1d;
    intNny = (dline.*ny)'*N1d;
    %HDG face matrices
    nodes1 = faceNodesApprox(iface,:);
     nodes2 = nodes1+nOfElementNodes; nodes3=nodes2+nOfElementNodes; nodes4=nodes3+nOfElementNodes;
    ind_face1 = (iface-1)*2*nOfFaceNodes + (1:nOfFaceNodes); 
     ind_face2 = ind_face1 + nOfFaceNodes;
    Auu(nodes1,nodes1) = Auu(nodes1,nodes1)+tau_f*Mf; 
     Auu(nodes2,nodes2) = Auu(nodes2,nodes2)+tau_f*Mf; 
    Aul(nodes1,ind_face1) = Aul(nodes1,ind_face1)-tau_f*Mf; 
     Aul(nodes2,ind_face2) = Aul(nodes2,ind_face2)-tau_f*Mf;
    ALl(nodes1,ind_face1) = ALl(nodes1,ind_face1)-Mfnx;
     ALl(nodes2,ind_face1) = ALl(nodes2,ind_face1)-Mfny;
     ALl(nodes3,ind_face2) = ALl(nodes3,ind_face2)-Mfnx;
     ALl(nodes4,ind_face2) = ALl(nodes4,ind_face2)-Mfny;    
    Apl(nodes1,ind_face1) = Apl(nodes1,ind_face1)-Mfnx; 
     Apl(nodes1,ind_face2) = Apl(nodes1,ind_face2)-Mfny;
    Arp(nodes1)=Arp(nodes1)+intN;
    All(ind_face1,ind_face1)=All(ind_face1,ind_face1)-tau_f*Mf;
     All(ind_face2,ind_face2)=All(ind_face2,ind_face2)-tau_f*Mf;
     %Arl([ind_face1,ind_face2])=Arl([ind_face1,ind_face2])+[intNnx,intNny];
 end


%%% Integrations over interface 

N1d=referenceElement.N1d;   %%%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  %%%%%%%%%%%%

%Initialization
ALutildaonI = zeros(4*nOfElementNodes,2*nOfFaceNodes);
AuuonI = zeros(2*nOfElementNodes,2*nOfElementNodes);
AuutildaonI = zeros(2*nOfElementNodes,2*nOfFaceNodes);
AputildaonI = zeros(nOfElementNodes,2*nOfFaceNodes);
Arutilda = zeros(1,2*nOfFaceNodes);
Autildautilda = zeros(2*nOfFaceNodes,2*nOfFaceNodes);
ge = zeros(2*nOfFaceNodes,1);



p=referenceElement.degree;
g=length(referenceElement.IPweights1d);
PtsInt=reshape(PtsInt',2,p+1)'; %interface p+1 nodes in the REFERENCE element
zg=referenceElement.N1d*PtsInt;%integration points on the interface in the REFERENCE element
shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zg);
Nzg = shapeFunctions(:,:,1)';  % Shape functions at integration points
shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,PtsInt);
NPtsInt = shapeFunctions(:,:,1)';  %1D Shape functions on the REFERENCE element at interface nodes
Pphy=NPtsInt*Xe; %p+1 interface nodes on the Physical element
Iprime=referenceElement.Nxi1d*Pphy;
I=referenceElement.N1d*Pphy; %integration pts on the interface


tau_f = tau(1);
dxdxi=Iprime(:,1);
dydxi=Iprime(:,2);

dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
dline = dxdxiNorm.*referenceElement.IPweights1d;
nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm;

%     figure(1); hold on;
%     quiver(I(:,1),I(:,2),nx,ny)

%Interface basic matrices
MfI = Nzg'*(spdiags(dline,0,ngf,ngf)*Nzg);         %for AuuonI
MfnxI = Nzg'*(spdiags(dline.*nx,0,ngf,ngf)*N1d);   %for AlutildaonI
MfnyI = Nzg'*(spdiags(dline.*ny,0,ngf,ngf)*N1d);   %for AlutildaonI
MffI = Nzg'*(spdiags(dline,0,ngf,ngf)*N1d);        %for AuutildaonI
MffII = N1d'*(spdiags(dline,0,ngf,ngf)*N1d);       %for Autildautilda
intNnxI = (dline.*nx)'*N1d;                        %for Arutilda
intNnyI = (dline.*ny)'*N1d;                        %for Arutilda
NnxIint = N1d'*[neumannu1(nx,ny,I,mu) .*dline];            %for ge  %curved interface wont work!!!!!
NnyIint = N1d'*[neumannu2(nx,ny,I,mu) .*dline];            %for ge

%XHDG interface matrices

ALutildaonI = ALutildaonI -[MfnxI    0*MfnxI ;
                MfnyI    MfnyI*0 ;
                MfnxI*0   MfnxI ;
                MfnyI*0   MfnyI ];  

AuuonI = [ MfI , MfI*0 ;
           MfI*0 , MfI ] .*tau_f;

AuutildaonI = AuutildaonI - [ MffI , MffI*0 ;
                              MffI*0 , MffI ] .*tau_f;  
            
AputildaonI = AputildaonI - [MfnxI , MfnyI] ; 

Autildautilda = Autildautilda - [MffII , MffII*0 ;
                                MffII*0 , MffII ].*tau_f;   

Arutilda = [intNnxI , intNnyI];

ge = [NnxIint ; NnyIint];

AutildaLonI=mu*ALutildaonI';
AutildaponI=-AputildaonI';   
AutildauonI=-AuutildaonI';

TL=[Autildautilda^-1]*[-AutildaLonI];
TU=[Autildautilda^-1]*[-AutildauonI];
TP=[Autildautilda^-1]*[-AutildaponI];
Ti=[Autildautilda^-1]*[ge];  

           
Alp=-Apl';
AlL=mu*ALl';
Alu=-Aul';


%%%%%%%%%%%%%%%%%%%%%%%%%
% % aux=ones(nOfElementNodes,1);
% % L_analy= [aux.*1 ; aux.*1; aux.*1; aux.*-1];
% % u_analy=analyticalVelocityStokes(Xe,mu);
% % aux2=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
% % gama=[analyticalVelocityStokes(Xe(aux2(1:nOfFaceNodes),:),mu);analyticalVelocityStokes(Xe(aux2(nOfFaceNodes+1:2*nOfFaceNodes),:),mu);...
% %    analyticalVelocityStokes(Xe(aux2(2*nOfFaceNodes+1:3*nOfFaceNodes),:),mu)];
% % p_analy=Xe(:,2)+Xe(:,1);
% % utilda_analy=analyticalVelocityStokes(Pphy,mu);
% % TestTest=ALL*L_analy+ALu*u_analy+ALl*gama+ALutildaonI*utilda_analy
% % Testa=AuL*L_analy+(Auu+AuuonI)*u_analy+Aup*p_analy
% % test=[AutildaLonI*L_analy+AutildaponI*p_analy+AutildauonI*u_analy+Autildautilda*utilda_analy-ge]
% % test45=utilda_analy-TL*L_analy-TP*p_analy-TU*u_analy+Ti;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ALL=ALL+ALutildaonI*TL;
ALu=ALu+ALutildaonI*TU;
ALp=ALutildaonI*TP;
AuL=AuL+AuutildaonI*TL;
Auu=Auu+AuuonI+AuutildaonI*TU;
Aup=Aup+AuutildaonI*TP;
ApL=AputildaonI*TL;
Apu=Apu+AputildaonI*TU;
App=AputildaonI*TP;



% Elemental mapping ( system A*[u L p lambda]=b - B*Lambda_e )
% lambda is the Lagrange multiplier to impose the constraint
% Lambda_e contains the trace of u at the faces AND the density mean (rho)
%A = zeros(7*nOfElementNodes+1);
A = zeros(7*nOfElementNodes);
indu = 1:2*nOfElementNodes; indL=2*nOfElementNodes + [1:4*nOfElementNodes]; 
indp = 6*nOfElementNodes + [1:nOfElementNodes];
A(indu,[indu,indL,indp]) = [Auu AuL Aup];
A(indL,[indu,indL,indp]) = [ALu ALL ALp];
A(indp,[indu,indL,indp]) = [Apu ApL App];
%A(end,indp) = Arp; A(indp,end) = Arp';
%b = zeros(7*nOfElementNodes+1,1);
b = zeros(7*nOfElementNodes,1);
b(indL) = -ALutildaonI*Ti;
b(indu) = fe-AuutildaonI*Ti;
b(indp) = -AputildaonI*Ti;
%Br = [zeros(7*nOfElementNodes,1) ; 1];
%B = [Aul ; ALl ; Apl ; zeros(1,6*nOfFaceNodes)];
B = [Aul ; ALl ; Apl]; 
%B = [Bl,Br]; %Lambda_e contains the face values of u and the rho value


% Br = [zeros(7*nOfElementNodes,1) ; 1];
% Bl = [Aul ; ALl ; Apl ; zeros(1,6*nOfFaceNodes)];
% B = [Bl,Br]; %Lambda_e contains the face values of u and the rho value
% 
%  LocalSolverMatElem=-A\B; LocalSolverMatElem=LocalSolverMatElem(1:end-1,:);
%  LocalSolverVecElem=A\b; LocalSolverVecElem=LocalSolverVecElem(1:end-1);

 
LocalSolverMatElem=-A\B; 
LocalSolverVecElem=A\b;  


% Updated Arl because of utilda on the interface; the new condition is
% [Arl][l] + [Arutilda][utilda]=0
%Arl=[Arl, 0]+[[Arutilda*TU], [Arutilda*TL], [Arutilda*TP]]*LocalSolverMatElem + [Arutilda*Ti]; 

%%%% Local Problem test %%%%% 
aux=zeros(nOfElementNodes,1);
%L_analy= [aux ;5.*Xe(:,2).^4; 5.*Xe(:,1).^4; aux];
L_analy= [aux.*1 ; aux.*1; aux.*1; aux.*-1];
u_analy=analyticalVelocityStokes(Xe,mu);
aux2=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
gama=[analyticalVelocityStokes(Xe(aux2(1:nOfFaceNodes),:),mu);analyticalVelocityStokes(Xe(aux2(nOfFaceNodes+1:2*nOfFaceNodes),:),mu);...
   analyticalVelocityStokes(Xe(aux2(2*nOfFaceNodes+1:3*nOfFaceNodes),:),mu)];
p_analy=Xe(:,2)+Xe(:,1);
utilda_analy=analyticalVelocityStokes(Pphy,mu);
% % % % rho=0;
% % % % % 
% % % % new_test=LocalSolverMatElem*[gama]+LocalSolverVecElem -[u_analy;L_analy;p_analy]
% % test1=ALL*L_analy+ALu*u_analy+ALp*p_analy+ALl*gama+ALutildaonI*Ti
% % test2=AuL*L_analy+Auu*u_analy+Aup*p_analy+Aul*gama-fe+AuutildaonI*Ti
% % test3=ApL*L_analy+Apu*u_analy+App*p_analy+Apl*gama+AputildaonI*Ti  



%ALL*L_analy+ALu*u_analy+ALl*gama+ALutildaonI*utilda_analy
%test=[AutildaLonI*L_analy+AutildaponI*p_analy+AutildauonI*u_analy+Autildautilda*utilda_analy-ge];




%save matrices A Aul Aql All fe sourceTerm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



