clear all; clc, %close all;
global R;

setpath

%mesh---------------------------------------------------------------------
degree=3;
meshName = ['mesh3_P',num2str(degree),'.mat']; 
if all(meshName(end-2:end)=='dcm')
    GenerateMatFileFromEZ4U(['Meshes/' meshName]);
end
load(['Meshes/' meshName(1:end-3) 'mat']);
X=2*X-1;

x=X(:,1); y=X(:,2);
% a = 0.51; b=0.03;
% LS = y-a+b*sin(4*pi*x);
centre=[0,0];
R = 0.45; 

LS = sqrt(x.^2+y.^2)-R;

figure(1),clf
plotMesh(X,T)
hold on, 
t = 0:0.01:2*pi;
[xx,yy]=meshgrid(t,t);
plot(R*cos(t),R*sin(t),'r-')
hold off
%---------------------------------------------------------------------------
referenceElement = createReferenceElementTri(degree);
Elements = SetElements(T,LS,[1,0],referenceElement);

%---------------------------------------------------------------------------
interfaceBasisAndQuadrature=computeBasisAndQuadratureAlongInterface(X,T,LS,referenceElement,Elements)
%Plot of nodes and normal vectors at integration points
hold on
for i=1:length(Elements.Int)
    quiver(interfaceBasisAndQuadrature.IPxy(:,1,i),interfaceBasisAndQuadrature.IPxy(:,2,i),interfaceBasisAndQuadrature.IPnormal(:,1,i),interfaceBasisAndQuadrature.IPnormal(:,2,i))
    plot(interfaceBasisAndQuadrature.nodesXY(:,1,i),interfaceBasisAndQuadrature.nodesXY(:,2,i),'o-')
end
hold off
%Computation of the length and the integral of f=x^2 on the circle
peri=0;
s1=0; s2=0;
fX=X(:,1).^2;
for i=1:length(Elements.Int)
   dl=interfaceBasisAndQuadrature.IPw(:,i); xy=interfaceBasisAndQuadrature.IPxy(:,:,i); 
   s1=s1+dl'*(xy(:,1).^2);
   peri=peri+sum(dl);
   Te=T(Elements.Int(i),:); fXe=fX(Te); fIP=interfaceBasisAndQuadrature.IPN(:,:,i)*fXe;
   s2=s2+dl'*fIP;
end
if(abs(peri-2*pi*R)>1.e-2) disp('Error in the length of the interface'), end
if(abs(s1-pi*R^3)>0.5e-2) disp('Error in the integral of f=x^2 on the interface (s1)'), end
if(abs(s2-pi*R^3)>0.5e-2) disp('Error in the integral of f=x^2 on the interface (s2)'), end
    

%---------------------------------------------------------------------------
AnalyticalSolution1=@(X) X(:,1).^2+X(:,2).^2;
AnalyticalSolution2=@(X) X(:,1)+X(:,2);

u1=AnalyticalSolution1(X); u2=AnalyticalSolution2(X); %CFE
nu1=1; nu2=2;

%___Plot of the jump of u on the interface (x,jump(u))
figure(11), clf, title('Plot of the jump of u at interface nodes')
%Numerical plot
hold on
for i=1:length(Elements.Int)
    Te=T(Elements.Int(i),:);
    N=interfaceBasisAndQuadrature.nodesN(:,:,i);
    u1I=N*u1(Te); u2I=N*u2(Te);
    jumpUe=u1I-u2I;
    xs=interfaceBasisAndQuadrature.nodesXY(:,1,i);
    ys=interfaceBasisAndQuadrature.nodesXY(:,2,i);
    theta=atan2(ys,xs); 
    %parche para +-pi
    k=find(abs(abs(theta)-pi)<1.e-10);if ~isempty(k)&any(abs(theta-theta(k))>pi/2), theta(k)=-theta(k); end
    plot(theta,jumpUe,'o-');
end
hold off
xlabel('theta'), ylabel('u1-u2')
%Analytical plot
theta=-pi:0.1:pi;
xs=R*cos(theta');ys=R*sin(theta');
jumpU=AnalyticalSolution1([xs,ys])-AnalyticalSolution2([xs,ys]);
hold on, plot(theta,jumpU,'r-'), hold off

%___Plot of the jump of the normal flux on the interface integration points
figure(12), clf, title('Plot of the jump of the normal flux at interface integration points')
%Numerical plot
hold on
for i=1:length(Elements.Int)
    Te=T(Elements.Int(i),:);
    dNdx=interfaceBasisAndQuadrature.IPdNdx(:,:,i);
    dNdy=interfaceBasisAndQuadrature.IPdNdy(:,:,i);
    jumpFlux=[dNdx*(nu1*u1(Te)-nu2*u2(Te)) dNdy*(nu1*u1(Te)-nu2*u2(Te))];
    nI=interfaceBasisAndQuadrature.IPnormal(:,:,i);
    jumpNormalFlux=nI(:,1).*jumpFlux(:,1)+nI(:,2).*jumpFlux(:,2);
    xs=interfaceBasisAndQuadrature.IPxy(:,1,i);
    ys=interfaceBasisAndQuadrature.IPxy(:,2,i);
    theta=atan2(ys,xs); 
    plot(theta,jumpNormalFlux,'-*');
end
hold off
xlabel('theta'), ylabel('jump(normal flux)')
%Analytical plot
theta=-pi:0.1:pi;
xs=R*cos(theta');ys=R*sin(theta');
jumpFlux=[nu1*2*xs-nu2,nu1*2*ys-nu2];
normal=-[xs/R,ys/R];
hold on, plot(theta,jumpFlux(:,1).*normal(:,1)+jumpFlux(:,2).*normal(:,2),'r-'), hold off

%___Plot of the jump of the normal flux on the interface NODES
figure(13), clf, title('Plot of the jump of the normal flux at interface NODES')
%Numerical plot
hold on
for i=1:length(Elements.Int)
    Te=T(Elements.Int(i),:);
    dNdx=interfaceBasisAndQuadrature.nodesdNdx(:,:,i);
    dNdy=interfaceBasisAndQuadrature.nodesdNdy(:,:,i);
    jumpFlux=[dNdx*(nu1*u1(Te)-nu2*u2(Te)) dNdy*(nu1*u1(Te)-nu2*u2(Te))];
    nI=interfaceBasisAndQuadrature.nodesNormal(:,:,i);
    jumpNormalFlux=nI(:,1).*jumpFlux(:,1)+nI(:,2).*jumpFlux(:,2);
    xs=interfaceBasisAndQuadrature.nodesXY(:,1,i);
    ys=interfaceBasisAndQuadrature.nodesXY(:,2,i);
    theta=atan2(ys,xs); 
    %parche para +-pi
    k=find(abs(abs(theta)-pi)<1.e-10);if ~isempty(k)&any(abs(theta-theta(k))>pi/2), theta(k)=-theta(k); end
    plot(theta,jumpNormalFlux,'-o');
end
hold off
xlabel('theta'), ylabel('jump(normal flux)')
%Analytical plot
theta=-pi:0.1:pi;
xs=R*cos(theta');ys=R*sin(theta');
jumpFlux=[nu1*2*xs-nu2,nu1*2*ys-nu2];
normal=-[xs/R,ys/R];
hold on, plot(theta,jumpFlux(:,1).*normal(:,1)+jumpFlux(:,2).*normal(:,2),'r-'), hold off

%___Plot of the normal fluxes on the interface NODES
figure(14), clf, title('Plot of the normal fluxes on interface NODES')
%Numerical plot
hold on
for i=1:length(Elements.Int)
    Te=T(Elements.Int(i),:);
    dNdx=interfaceBasisAndQuadrature.nodesdNdx(:,:,i);
    dNdy=interfaceBasisAndQuadrature.nodesdNdy(:,:,i);
    Flux1=[dNdx*(nu1*u1(Te)) dNdy*(nu1*u1(Te))];
    Flux2=[dNdx*(nu2*u2(Te)) dNdy*(nu2*u2(Te))];
    nI=interfaceBasisAndQuadrature.nodesNormal(:,:,i);
    normalFlux1=nI(:,1).*Flux1(:,1)+nI(:,2).*Flux1(:,2);
    normalFlux2=nI(:,1).*Flux2(:,1)+nI(:,2).*Flux2(:,2);
    xs=interfaceBasisAndQuadrature.nodesXY(:,1,i);
    ys=interfaceBasisAndQuadrature.nodesXY(:,2,i);
    theta=atan2(ys,xs); 
    %parche para +-pi
    k=find(abs(abs(theta)-pi)<1.e-10);if ~isempty(k)&any(abs(theta-theta(k))>pi/2), theta(k)=-theta(k); end
    plot(theta,[normalFlux1,normalFlux2],'-o');
end
legend('normal flux 1','-normal flux 2')
hold off
xlabel('theta'), ylabel('normal flux')
%Analytical plot
theta=-pi:0.1:pi;
xs=R*cos(theta');ys=R*sin(theta');
jumpFlux=[nu1*2*xs-nu2,nu1*2*ys-nu2];
normal=-[xs/R,ys/R];
hold on, plot(theta,[(nu1*2*xs).*normal(:,1)+(nu1*2*ys).*normal(:,2),nu2*(normal(:,1)+normal(:,2))],'r-'), hold off

%___Integral of (u1-u2)^2 and jump(normal flux)^2
integralJumpU2=0; integralJumpNormalFlux2=0;
for i=1:length(Elements.Int)
    Te=T(Elements.Int(i),:);
    dNdx=interfaceBasisAndQuadrature.IPdNdx(:,:,i);
    dNdy=interfaceBasisAndQuadrature.IPdNdy(:,:,i);
    jumpFlux=[dNdx*(nu1*u1(Te)-nu2*u2(Te)) dNdy*(nu1*u1(Te)-nu2*u2(Te))];
    nI=interfaceBasisAndQuadrature.IPnormal(:,:,i);
    jumpNormalFlux=nI(:,1).*jumpFlux(:,1)+nI(:,2).*jumpFlux(:,2);
    dl=interfaceBasisAndQuadrature.IPw(:,i);
    integralJumpNormalFlux2=integralJumpNormalFlux2+dl'*(jumpNormalFlux.^2);
    N=interfaceBasisAndQuadrature.IPN(:,:,i);
    u1I=N*u1(Te); u2I=N*u2(Te);
    jumpUe=u1I-u2I;
    integralJumpU2 = integralJumpU2 + dl'*(jumpUe.^2);
end
normL2jumpU=sqrt(integralJumpU2)
normL2jumpU_analytical=0.8297576161
normL2jumpNormalFlux=sqrt(integralJumpNormalFlux2)
normL2jumpNormalFlux_analytical=3.687811628
