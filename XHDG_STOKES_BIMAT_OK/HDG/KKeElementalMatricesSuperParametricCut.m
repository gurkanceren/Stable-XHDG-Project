function [LocalSolverMatElem,LocalSolverVecElem,Alu,AlL,Alp,All,Arl,AI] = KKeElementalMatricesSuperParametricCut(Xe,referenceElement,tau,LSe)
global mu1 mu2;
nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
faceNodesApprox = referenceElement.faceNodes;
faceNodes = referenceElement.faceNodes;
nOfFaces = 3; %triangles


[Nold,Nxi_old,Neta_old,N,Nxi,Neta,weight,ngauss,mu,PtsInt,FaceInfo]=modifiedReferenceElement(referenceElement,LSe);
referenceElement.IPweights=weight';

J11 = Nxi_old*Xe(:,1); J12 = Nxi_old*Xe(:,2); 
J21 = Neta_old*Xe(:,1); J22 = Neta_old*Xe(:,2); 
detJ = J11.*J22-J12.*J21;

%maybe we should use bsxfun instead of diagonal matrices...
%dvolu = spdiags(referenceElement.IPweights.*detJ(1:ngauss),0,ngauss,ngauss);  % !!!
dvolu = spdiags(weight'.*detJ,0,ngauss,ngauss);  % !!!
dvoluaux = spdiags(weight'.*detJ.*-mu,0,ngauss,ngauss);  % !!!

invJ11 = spdiags(J22./detJ,0,ngauss,ngauss);
invJ12 = spdiags(-J12./detJ,0,ngauss,ngauss);
invJ21 = spdiags(-J21./detJ,0,ngauss,ngauss);
invJ22 = spdiags(J11./detJ,0,ngauss,ngauss);
% xy-derivatives for approximation
Nx = invJ11*Nxi + invJ12*Neta;
Ny = invJ21*Nxi + invJ22*Neta;

%Computation of r.h.s. source term (analytical laplacian)
Xg = Nold*Xe;
% figure(1); hold on;
% plot(Xg(1:(ngauss/2),1),Xg(1:(ngauss/2),2),'*b')
sourceTerm = sourceTermStokes(Xg); %Missing!!!
fe = N'*(dvolu*sourceTerm);
fe = reshape(fe,4*nOfElementNodes,1);
%Basic elemental matrices 
Me = N'*(dvolu*N); Cxe=N'*(dvolu*Nx); Cye=N'*(dvolu*Ny);  
Cxeaux=Nx'*(dvoluaux*N); Cyeaux=Ny'*(dvoluaux*N);
%HDG elemental matrices
aux1=1:2*nOfElementNodes; aux2=aux1+2*nOfElementNodes; aux3=aux2+2*nOfElementNodes; aux4=aux3+2*nOfElementNodes;
ALL = zeros(8*nOfElementNodes,8*nOfElementNodes);
 ALL(aux1,aux1)=Me; ALL(aux2,aux2)=Me; ALL(aux3,aux3)=Me; ALL(aux4,aux4)=Me; 
ALu = zeros(8*nOfElementNodes,4*nOfElementNodes);
 ALu(aux1,aux1)=Cxe'; ALu(aux2,aux1)=Cye'; ALu(aux3,aux2)=Cxe'; ALu(aux4,aux2)=Cye';
AuL = zeros(4*nOfElementNodes,8*nOfElementNodes); 
AuL(aux1,aux1)=Cxeaux'; AuL(aux1,aux2)=Cyeaux'; AuL(aux2,aux3)=Cxeaux'; AuL(aux2,aux4)=Cyeaux';
Aup = [Cxe;Cye]; Apu = Aup';


%1D Numerical quadrature
N1d = referenceElement.N1d; Nx1d = referenceElement.Nxi1d;
IPw_f = referenceElement.IPweights1d; ngf = length(IPw_f);

%% Faces computations
ALl = zeros(8*nOfElementNodes,3*4*nOfFaceNodes);
AlL = zeros(3*4*nOfFaceNodes,8*nOfElementNodes);
Auu = zeros(4*nOfElementNodes,4*nOfElementNodes);
Aul = zeros(4*nOfElementNodes,3*4*nOfFaceNodes);
Apl = zeros(2*nOfElementNodes,3*4*nOfFaceNodes);
All = zeros(3*4*nOfFaceNodes,3*4*nOfFaceNodes);
Arp = zeros(1,2*nOfElementNodes);
Arl = zeros(1,3*4*nOfFaceNodes);


 for iface = 1:nOfFaces
   
    if ~isempty(FaceInfo)
        k=find(FaceInfo(:,1)==iface);
    else
        k=[];
    end
    
    Xf = Xe(faceNodes(iface,:),:); % Nodes in the face
    nodes = faceNodes(iface,:);
    ind_face =(iface-1)*(4*nOfFaceNodes) + (1:(4*nOfFaceNodes));
    if ~isempty(k);  %Cut face 
        
        %Shape Functions for cut faces:
        [N_old,Nx_cutf_old,Nx_cutf,N_cutf,wgp_f,ngf]=modifiedReferenceElement1Cut(referenceElement,LSe,nodes);
        %Integration weight
        Xgf=N_old*Xf;
        N1denr=N_cutf;
        N1d=N_cutf;
        IPw_f=wgp_f;
        dxdxi = Nx_cutf_old*Xf(:,1); dydxi = Nx_cutf_old*Xf(:,2);
        mul=ones((ngf/2),1)*mu1;
        mur=ones((ngf/2),1)*mu2;
        mu=[mul;mur];
        %gama(ind_face) = [analyticalVelocityStokesH(Xf)]; %for TEST
        
    elseif isempty(k) && any(LSe(referenceElement.faceNodes(iface,:))>0) ;  %Material face 1
        N1d_old = referenceElement.N1d;    Nx1d_old = referenceElement.Nxi1d;
        Xgf= N1d_old*Xf;
%         figure(1); hold on;
%         plot(Xgf(:,1),Xgf(:,2),'md')
        N1d= [N1d_old  N1d_old*0] ;
        N1denr= [N1d_old  N1d_old*1];
        IPw_f = referenceElement.IPweights1d; ngf = length(IPw_f);
        dxdxi = Nx1d_old*Xf(:,1); dydxi = Nx1d_old*Xf(:,2);
        mu=ones((ngf),1)*mu1;
        %gama(ind_face) = [analyticalVelocityStokes(Xf)]; %for TEST
        %gama(ind_face) = [Xf(:,1)'.*0+10,Xf(:,1)'.*0,Xf(:,1)'.*0+10,Xf(:,1)'.*0]; %for TEST
        
    elseif  isempty(k) && any(LSe(referenceElement.faceNodes(iface,:))<0) ;  %Material face 2
        
        N1d_old = referenceElement.N1d;    Nx1d_old = referenceElement.Nxi1d;
        Xgf= N1d_old*Xf;
%         figure(1); hold on;
%         plot(Xgf(:,1),Xgf(:,2),'cd')
        N1d= [N1d_old  N1d_old*0] ;
        N1denr= [N1d_old  N1d_old*-1];
        IPw_f = referenceElement.IPweights1d; ngf = length(IPw_f);
        dxdxi = Nx1d_old*Xf(:,1); dydxi = Nx1d_old*Xf(:,2);
        mu=ones((ngf),1)*mu2;
        %gama(ind_face) = [analyticalVelocityStokes(Xf)]; %for TEST
        %gama(ind_face) = [Xf(:,1)'.*0+4,Xf(:,1)'.*0,Xf(:,1)'.*0+4,Xf(:,1)'.*0]; %for TEST
    end
   
    tau_f = tau(iface);
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = dxdxiNorm.*IPw_f;
    nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm;

    %Face basic matrices
    Mf = N1denr'*(spdiags(mu.*dline,0,ngf,ngf)*N1denr);
    Mff = N1denr'*(spdiags(mu.*dline,0,ngf,ngf)*N1d);
    Mfff = N1d'*(spdiags(mu.*dline,0,ngf,ngf)*N1d);
    
    Mfnx = N1d'*(spdiags(dline.*nx,0,ngf,ngf)*N1d);
    Mffnx = N1denr'*(spdiags(dline.*nx,0,ngf,ngf)*N1d);
    Mffnxaux = N1d'*(spdiags(-mu.*dline.*nx,0,ngf,ngf)*N1denr);
    
    Mfny = N1d'*(spdiags(dline.*ny,0,ngf,ngf)*N1d);
    Mffny = N1denr'*(spdiags(dline.*ny,0,ngf,ngf)*N1d);
    Mffnyaux = N1d'*(spdiags(-mu.*dline.*ny,0,ngf,ngf)*N1denr);
    
    intN = dline'*N1denr;
    intNnx = (dline.*nx)'*N1d;
    intNny = (dline.*ny)'*N1d;
    %HDG face matrices
    nodes1 = faceNodesApprox(iface,:);
    nodes2 = nodes1+nOfElementNodes; 
    nodes3=  nodes2+nOfElementNodes;
    nodes4=  nodes3+nOfElementNodes;
    nodes5=  nodes4+nOfElementNodes;
    nodes6=  nodes5+nOfElementNodes;
    nodes7=  nodes6+nOfElementNodes;

    nodes1aux = [nodes1 nodes1+nOfElementNodes];
    nodes2aux = [nodes3 nodes3+nOfElementNodes];
    nodes3aux = [nodes5 nodes5+nOfElementNodes];
    nodes4aux = [nodes7 nodes7+nOfElementNodes];
    nodes1u=[nodes1 nodes1+nOfElementNodes];
    nodes2u=[nodes3 nodes3+nOfElementNodes];
    nodes1p=[nodes1 nodes1+nOfElementNodes];
    
    ind_face1 = (iface-1)*4*nOfFaceNodes + (1:2*nOfFaceNodes);
    ind_face2 = ind_face1 + 2*nOfFaceNodes;
    Auu(nodes1u,nodes1u) = Auu(nodes1u,nodes1u)+tau_f*Mf;
    Auu(nodes2u,nodes2u) = Auu(nodes2u,nodes2u)+tau_f*Mf;
    
    Aul(nodes1u,ind_face1) = Aul(nodes1u,ind_face1)-tau_f*Mff;
    Aul(nodes2u,ind_face2) = Aul(nodes2u,ind_face2)-tau_f*Mff;
    
    ALl(nodes1aux,ind_face1) = ALl(nodes1aux,ind_face1)-Mffnx;
    ALl(nodes2aux,ind_face1) = ALl(nodes2aux,ind_face1)-Mffny;
    ALl(nodes3aux,ind_face2) = ALl(nodes3aux,ind_face2)-Mffnx;
    ALl(nodes4aux,ind_face2) = ALl(nodes4aux,ind_face2)-Mffny;
    
    AlL(ind_face1,nodes1aux) =  AlL(ind_face1,nodes1aux)+Mffnxaux;
    AlL(ind_face1,nodes2aux) = AlL(ind_face1,nodes2aux)+Mffnyaux;
    AlL(ind_face2,nodes3aux) = AlL(ind_face2,nodes3aux)+Mffnxaux;
    AlL(ind_face2,nodes4aux) = AlL(ind_face2,nodes4aux)+Mffnyaux;
    
    Apl(nodes1p,ind_face1) = Apl(nodes1p,ind_face1)-Mffnx;
    Apl(nodes1p,ind_face2) = Apl(nodes1p,ind_face2)-Mffny;
    Arp(nodes1p)=Arp(nodes1p)+intN;
    All(ind_face1,ind_face1)=All(ind_face1,ind_face1)-tau_f*Mfff;
    All(ind_face2,ind_face2)=All(ind_face2,ind_face2)-tau_f*Mfff;
    Arl([ind_face1,ind_face2])=Arl([ind_face1,ind_face2])+[intNnx,intNny];
 end

Alp=-Apl';
Alu=-Aul';

%%% Integrations over interface 

N1d=referenceElement.N1d;   %%%%%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  %%%%%%%%%%%%

%Initialization
ALutildaonI = zeros(8*nOfElementNodes,2*nOfFaceNodes);
AutildaLonI = zeros(2*nOfFaceNodes,8*nOfElementNodes);
AuuonI = zeros(4*nOfElementNodes,4*nOfElementNodes);
AuutildaonI = zeros(4*nOfElementNodes,2*nOfFaceNodes);
AputildaonI = zeros(2*nOfElementNodes,2*nOfFaceNodes);
Arutilda = zeros(1,2*nOfFaceNodes);
Autildautilda = zeros(2*nOfFaceNodes,2*nOfFaceNodes);

newterm2= zeros(8*nOfElementNodes,1);
newterm2_y1 = zeros(2*nOfElementNodes,1);
newterm2_y2 = zeros(2*nOfElementNodes,1);
newterm2_x1 = zeros(2*nOfElementNodes,1);
newterm2_x2 = zeros(2*nOfElementNodes,1);
newterm3 = zeros(2*nOfElementNodes,1);

newterm11=zeros(2*nOfElementNodes,1);
newterm12=zeros(2*nOfElementNodes,1);
newterm1=zeros(4*nOfElementNodes,1);
AI=0;
ge1=zeros(nOfFaceNodes,1);
ge2=zeros(nOfFaceNodes,1);
ge=zeros(2*nOfFaceNodes,1);

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
ngf=size(referenceElement.IPweights1d,1);

    NzgR=[Nzg Nzg*-1];
    NzgL=[Nzg Nzg*1];
    
    muL=ones(g,1)*mu1;
    muR=ones(g,1)*mu2;
   
for igauss=1:g 
    
         [sD1 , sD2]=calculateSD(I(igauss,:));
         newterm11 = newterm11+ tau_f*(NzgR(igauss,:)'*-sD1 *muR(igauss))*dline(igauss);
         newterm12 = newterm12+ tau_f*(NzgR(igauss,:)'*-sD2 *muR(igauss))*dline(igauss);
         newterm2_y1 = newterm2_y1 - NzgR(igauss,:)'*-ny(igauss)*sD1*dline(igauss);
         newterm2_y2 = newterm2_y2 - NzgR(igauss,:)'*-ny(igauss)*sD2*dline(igauss);
         newterm2_x1 = newterm2_x1 - NzgR(igauss,:)'*-nx(igauss)*sD1*dline(igauss);
         newterm2_x2 = newterm2_x2 - NzgR(igauss,:)'*-nx(igauss)*sD2*dline(igauss); 
         newterm3 = newterm3 - NzgR(igauss,:)'*[sD1 sD2]*[-nx(igauss); -ny(igauss)]*dline(igauss);
         AI=AI+[sD1 sD2]*[-nx(igauss); -ny(igauss)]*dline(igauss);
         ge1 = ge1 + [tau_f*N1d(igauss,:)'*(muR(igauss)*sD1)*dline(igauss)]+ [N1d(igauss,:)'*NeumannCondition1([nx(igauss) ny(igauss)],[I(igauss,1),I(igauss,2)])*dline(igauss)];
         ge2 = ge2 + [tau_f*N1d(igauss,:)'*(muR(igauss)*sD2)*dline(igauss)]+ [N1d(igauss,:)'*NeumannCondition2([nx(igauss) ny(igauss)],[I(igauss,1),I(igauss,2)])*dline(igauss)]; 
         
end
    
newterm1=[newterm11;newterm12];
newterm2 = [newterm2_x1;newterm2_y1;newterm2_x2;newterm2_y2];
ge=[ge1;ge2]; 
%Interface basic matrices
MfI = NzgR'*(spdiags(muR.*dline,0,ngf,ngf)*NzgR)+NzgL'*(spdiags(muL.*dline,0,ngf,ngf)*NzgL);    %for AuuonI
MfnxI = (NzgL-NzgR)'*(spdiags(dline.*nx,0,ngf,ngf)*N1d);  %-(NzgR)'*(spdiags(dline.*nx,0,ngf,ngf)*N1d) ;   %for AlutildaonI
MfnyI = (NzgL-NzgR)'*(spdiags(dline.*ny,0,ngf,ngf)*N1d);  %-(NzgR)'*(spdiags(dline.*ny,0,ngf,ngf)*N1d);   %for AlutildaonI

MfnxIaux = N1d'*(spdiags(muL.*dline.*nx,0,ngf,ngf)*NzgL)-N1d'*(spdiags(muR.*dline.*nx,0,ngf,ngf)*NzgR);
MfnyIaux = N1d'*(spdiags(muL.*dline.*ny,0,ngf,ngf)*NzgL)-N1d'*(spdiags(muR.*dline.*ny,0,ngf,ngf)*NzgR);

MffI =NzgR'*(spdiags(muR.*dline,0,ngf,ngf)*N1d)+NzgL'*(spdiags(muL.*dline,0,ngf,ngf)*N1d);        %for AuutildaonI

MffII = N1d'*(spdiags((muL+muR).*dline,0,ngf,ngf)*N1d);       %for Autildautilda
intNnxI = (dline.*nx)'*N1d;                        %for Arutilda
intNnyI = (dline.*ny)'*N1d;                        %for Arutilda

%XHDG interface matrices

ALutildaonI = ALutildaonI -[MfnxI    0*MfnxI ;
    MfnyI    MfnyI*0 ;
    MfnxI*0   MfnxI ;
    MfnyI*0   MfnyI ];
            
AutildaLonI = AutildaLonI -[MfnxIaux    MfnyIaux   MfnyIaux.*0  MfnyIaux.*0;
                             MfnxIaux.*0    MfnyIaux.*0   MfnxIaux  MfnyIaux];
             

AuuonI = [ MfI , MfI*0 ;
           MfI*0 , MfI ] .*tau_f;

AuutildaonI = AuutildaonI - [ MffI , MffI*0 ;
                              MffI*0 , MffI ] .*tau_f;  
            
AputildaonI = AputildaonI - [MfnxI , MfnyI] ; 

Autildautilda = Autildautilda - [MffII , MffII*0 ;
                                MffII*0 , MffII ].*tau_f;   



AutildaponI=-AputildaonI';   
AutildauonI=-AuutildaonI';

TL=[Autildautilda^-1]*[-AutildaLonI];
TU=[Autildautilda^-1]*[-AutildauonI];
TP=[Autildautilda^-1]*[-AutildaponI];
be=[Autildautilda^-1]*ge;          

    
    Cu=AuutildaonI*be;
    Cq=ALutildaonI*be;
    Cp=AputildaonI*be;


ALL=ALL+ALutildaonI*TL;
ALu=ALu+ALutildaonI*TU;
ALp=ALutildaonI*TP;
AuL=AuL+AuutildaonI*TL;
Auu=Auu+AuuonI+AuutildaonI*TU;
Aup=Aup+AuutildaonI*TP;
ApL=AputildaonI*TL;
Apu=Apu+AputildaonI*TU;
App=AputildaonI*TP;

%p_analy=analyticalP(Xe,mu1,mu2);
% p_analy=[p_analy;p_analy.*0];
% Elemental mapping ( system A*[u L p lambda]=b - B*Lambda_e )
% lambda is the Lagrange multiplier to impose the constraint
% Lambda_e contains the trace of u at the faces AND the density mean (rho)
A = zeros(14*nOfElementNodes+1);
indu = 1:4*nOfElementNodes; indL=4*nOfElementNodes + [1:8*nOfElementNodes]; 
indp = 12*nOfElementNodes + [1:2*nOfElementNodes];
A(indu,[indu,indL,indp]) = [Auu AuL Aup];
A(indL,[indu,indL,indp]) = [ALu ALL ALp];
A(indp,[indu,indL,indp]) = [Apu ApL App];
A(end,indp) = Arp; A(indp,end) = Arp';
b = zeros(14*nOfElementNodes+1,1);
b(indL) = zeros(8*nOfElementNodes,1)-newterm2-Cq;
b(indu) = fe-newterm1-Cu;
b(indp) = zeros(2*nOfElementNodes,1)-newterm3-Cp;
Br = [zeros(14*nOfElementNodes,1); 1];
Bl = [Aul ; ALl ; Apl ; zeros(1,12*nOfFaceNodes)];
B=[Bl,Br];


LocalSolverMatElem=-A\B; LocalSolverMatElem=LocalSolverMatElem(1:end-1,:);
LocalSolverVecElem=A\b;  LocalSolverVecElem=LocalSolverVecElem(1:end-1);


% L_analy=analyticalLStokesH(Xe);
% L_analy=reshape(L_analy,nOfElementNodes*8,1);
% u_analy=analyticalVelocityStokesH(Xe);
% u_analy=reshape(u_analy,nOfElementNodes*4,1);
% p_analy=analyticalP(Xe);
% utilda_analy=analyticalVelocityStokes(Pphy);
% utilda_analy=reshape(utilda_analy,nOfElementNodes*2,1);

% 
% new_test=LocalSolverMatElem*[gama 0]'+LocalSolverVecElem -[u_analy;L_analy;p_analy];
% test1=ALL*L_analy+ALu*u_analy+ALl*gama'+ ALp*p_analy 
% test2=AuL*L_analy+Auu*u_analy+Aup*p_analy+Aul*gama'-fe
% test3=ApL*L_analy+Apu*u_analy+App*p_analy+Apl*gama'

% Updated Arl because of utilda on the interface; the new condition is
% [Arl][l] + [Arutilda][utilda]=0
%Arl=[Arl, 0]+[[Arutilda*TU], [Arutilda*TL], [Arutilda*TP]]*LocalSolverMatElem + [Arutilda*Ti]; 


%%% Local Problem test %%%%% 

 %%x=Xe(:,1); y=Xe(:,2);
% % re1=1/mu1;
% % re2=1/mu2;
% % Lmbd1=(re1/2)-sqrt(((re1^2)/4)+4*pi^2);
% % Lmbd2=(re2/2)-sqrt(((re2^2)/4)+4*pi^2);
% % aux=y.*0;
% % 
% % %p_analy=[x+y; (x+y).*0];
% % 
% % for i=1:length(y)
% %     
% %     if y<= 0
% % mu=mu1;
% %     else
% % mu=mu2;
% %     end
% %  
% % L_analy=[(2.*x.*y)/mu ; zeros(nOfElementNodes,1); (x.^2)/mu ; ...
% %     zeros(nOfElementNodes,1); (-y.^2)/mu ;zeros(nOfElementNodes,1); (-2.*x.*y)/mu ;zeros(nOfElementNodes,1)];
% % 
% % u_analy=analyticalVelocityStokesH(Xe);
% % u_analy=reshape(u_analy,nOfElementNodes*4,1);
% % %p_analy=analyticalP(Xe);
% % p_analy=[x+y; (x+y).*0];
% % 
% % ALL*L_analy+ALu*u_analy+ALl*gama'+ALp*p_analy;
% % AuL*L_analy+Auu*u_analy+Aup*p_analy+Aul*gama'-fe;
% % ApL*L_analy+Apu*u_analy+App*p_analy+Apl*gama';
% % 
% %test=[-AutildaLonI*L_analy+AutildaponI*p_analy+AutildauonI*u_analy+Autildautilda*utilda_analy-ge];
% test1=ALL*L_analy+ALu*u_analy+ALl*gama'
% test2=AuL*L_analy+Auu*u_analy+Aup*p_analy+Aul*gama'-fe
% test3=ApL*p_analy+Apu*u_analy+App*p_analy+Apl*gama'  %%% EQN (18) in the paper
% 
% % %save matrices A Aul Aql All fe sourceTerm
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



