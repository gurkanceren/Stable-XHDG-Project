function [LocalSolverMatElem,LocalSolverVecElem,Alu,AlL,Alp,All,Arl] = KKeElementalMatricesSuperParametric(mu,Xe,referenceElement,tau)


nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
faceNodesApprox = referenceElement.faceNodes;
faceNodesGeo = referenceElement.faceNodesGeo;
nOfFaces = 3; %triangles

% Information of the reference element
N = referenceElement.N;
Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
NxiGeo = referenceElement.NxiGeo; NetaGeo = referenceElement.NetaGeo;
NGeo = referenceElement.NGeo;
N1d = referenceElement.N1d; %Nx1d = referenceElement.Nxi1d;
Nx1dGeo = referenceElement.N1dxiGeo;
%Numerical quadrature
IPw_f = referenceElement.IPweights1d; ngf = length(IPw_f);
IPw = referenceElement.IPweights; ngauss = length(IPw);

%Void Elements set to zero
%if ~isempty(d2) IPw_f=IPw_f.*0;  IPw=IPw.*0; end
    
%%Volume computations
% Jacobian
J11 = NxiGeo*Xe(:,1); J12 = NxiGeo*Xe(:,2); 
J21 = NetaGeo*Xe(:,1); J22 = NetaGeo*Xe(:,2); 
detJ = J11.*J22-J12.*J21;
%maybe we should use bsxfun instead of diagonal matrices...
dvolu = spdiags(referenceElement.IPweights.*detJ,0,ngauss,ngauss);
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
% plot(Xg(:,1),Xg(:,2),'*y')
sourceTerm = sourceTermStokes(Xg); %Missing!!!
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

%% Faces computations
ALl = zeros(4*nOfElementNodes,3*2*nOfFaceNodes);
Auu = zeros(2*nOfElementNodes,2*nOfElementNodes);
Aul = zeros(2*nOfElementNodes,3*2*nOfFaceNodes);
Apl = zeros(nOfElementNodes,3*2*nOfFaceNodes);
All = zeros(3*2*nOfFaceNodes,3*2*nOfFaceNodes);
Arp = zeros(1,nOfElementNodes);
Arl = zeros(1,3*2*nOfFaceNodes);
for iface = 1:nOfFaces    
    taumu_f = tau(iface)*mu;
    Xf = Xe(faceNodesGeo(iface,:),:); % Nodes in the face
    dxdxi = Nx1dGeo*Xf(:,1); dydxi = Nx1dGeo*Xf(:,2);
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = dxdxiNorm.*IPw_f;
    nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm;
    %Face basic matrices
    Mf = N1d'*(spdiags(dline,0,ngf,ngf)*N1d);
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
    Auu(nodes1,nodes1) = Auu(nodes1,nodes1)+taumu_f*Mf; 
     Auu(nodes2,nodes2) = Auu(nodes2,nodes2)+taumu_f*Mf; 
    Aul(nodes1,ind_face1) = Aul(nodes1,ind_face1)-taumu_f*Mf; 
     Aul(nodes2,ind_face2) = Aul(nodes2,ind_face2)-taumu_f*Mf;
    ALl(nodes1,ind_face1) = ALl(nodes1,ind_face1)-Mfnx;
     ALl(nodes2,ind_face1) = ALl(nodes2,ind_face1)-Mfny;
     ALl(nodes3,ind_face2) = ALl(nodes3,ind_face2)-Mfnx;
     ALl(nodes4,ind_face2) = ALl(nodes4,ind_face2)-Mfny;    
    Apl(nodes1,ind_face1) = Apl(nodes1,ind_face1)-Mfnx; 
     Apl(nodes1,ind_face2) = Apl(nodes1,ind_face2)-Mfny;
    Arp(nodes1)=Arp(nodes1)+intN;
    All(ind_face1,ind_face1)=All(ind_face1,ind_face1)-taumu_f*Mf;
     All(ind_face2,ind_face2)=All(ind_face2,ind_face2)-taumu_f*Mf;
    Arl([ind_face1,ind_face2])=Arl([ind_face1,ind_face2])+[intNnx,intNny];
end
Alp=-Apl';
AlL=mu*ALl';
Alu=-Aul';
% p_analy=analyticalP(Xe,mu,mu);
% p_analy=p_analy(1:nOfElementNodes);
% Elemental mapping ( system A*[u L p lambda]=b - B*Lambda_e )
% lambda is the Lagrange multiplier to impose the constraint
% Lambda_e contains the trace of u at the faces AND the density mean (rho)
A = zeros(7*nOfElementNodes+1);
indu = 1:2*nOfElementNodes; indL=2*nOfElementNodes + (1:4*nOfElementNodes); 
indp = 6*nOfElementNodes + (1:nOfElementNodes);
A(indu,[indu,indL,indp]) = [Auu AuL Aup];
A(indL,[indu,indL]) = [ALu ALL];
A(indp,indu) = Apu;
A(end,indp) = Arp; A(indp,end) = Arp';
b = zeros(7*nOfElementNodes+1,1);
b(indu) = fe;
Br = [zeros(7*nOfElementNodes,1) ;  1];
Bl = [Aul ; ALl ; Apl ; zeros(1,6*nOfFaceNodes)];
B = [Bl,Br]; %Lambda_e contains the face values of u and the rho value

LocalSolverMatElem=-A\B; LocalSolverMatElem=LocalSolverMatElem(1:end-1,:);
LocalSolverVecElem=A\b; LocalSolverVecElem=LocalSolverVecElem(1:end-1);


% % %%%%% Local problem test  %%%%%% 
% aux=ones(nOfElementNodes,1);
% L_analy= [aux; aux ; aux; aux.*-1];
% u_analy=analyticalVelocityStokes(Xe,mu);
% % aux2=reshape(referenceElement.faceNodes',size(referenceElement.faceNodes,1)*size(referenceElement.faceNodes,2),1);
% % gama=[analyticalVelocityStokesH(Xe(aux2(1:nOfFaceNodes),:),mu);analyticalVelocityStokesH(Xe(aux2(nOfFaceNodes+1:2*nOfFaceNodes),:),mu);...
% %    analyticalVelocityStokesH(Xe(aux2(2*nOfFaceNodes+1:3*nOfFaceNodes),:),mu)];
% p_analy=Xe(:,2)+Xe(:,1);
% 
% rho=Arp*p_analy;
% 
% % x=Xe(:,1); y=Xe(:,2);
% % L_analy=[(2.*x.*y)/mu ; (x.^2)/mu ; (-y.^2)/mu;(-2.*x.*y)/mu];
% % u_analy=analyticalVelocityStokesH(Xe);
% % u_analy=reshape(u_analy,nOfElementNodes*4,1);
% % p_analy=x+y;
% % u_analy=[u_analy(1:6); u_analy(13:18)];
% % test1=ALL*L_analy+ALu*u_analy+ALl*gama;
% % test2=Apu*u_analy+Apl*gama;
% % test3=Auu*u_analy+AuL*L_analy+Aup*p_analy+Aul*gama-fe;

%save matrices A Aul Aql All fe sourceTerm



