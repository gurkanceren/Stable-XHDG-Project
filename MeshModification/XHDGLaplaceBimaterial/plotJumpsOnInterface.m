%___Plot of the jump of u on the interface (x,jump(u))
figure(11), clf, title('Plot of the jump of u at interface nodes')
%Numerical plot
hold on
for i=1:length(Elements.Int)
    ielem=Elements.Int(i);
    Te=(ielem-1)*nOfElementNodes + (1:nOfElementNodes);
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

%___Plot of the jump of the normal flux on the interface integration points
figure(12), clf, title('Plot of the jump of the normal flux at interface integration points')
%Numerical plot
hold on
for i=1:length(Elements.Int)
    ielem=Elements.Int(i);
    Te=(ielem-1)*nOfElementNodes + (1:nOfElementNodes);
    N=interfaceBasisAndQuadrature.IPN(:,:,i);
    nI=interfaceBasisAndQuadrature.IPnormal(:,:,i);
    jumpFlux=[N*(qx1(Te)-qx2(Te)) N*(qy1(Te)-qy2(Te))];
    jumpNormalFlux=nI(:,1).*jumpFlux(:,1)+nI(:,2).*jumpFlux(:,2);
    xs=interfaceBasisAndQuadrature.IPxy(:,1,i);
    ys=interfaceBasisAndQuadrature.IPxy(:,2,i);
    theta=atan2(ys,xs); 
    plot(theta,jumpNormalFlux,'-*');
end
hold off
xlabel('theta'), ylabel('jump(normal flux)')

%___Plot of the jump of the normal flux on the interface NODES
figure(13), clf, title('Plot of the jump of the normal flux at interface NODES')
%Numerical plot
hold on
for i=1:length(Elements.Int)
    ielem=Elements.Int(i);
    Te=(ielem-1)*nOfElementNodes + (1:nOfElementNodes);
    N=interfaceBasisAndQuadrature.nodesN(:,:,i);
    nI=interfaceBasisAndQuadrature.nodesNormal(:,:,i);
    jumpFlux=[N*(qx1(Te)-qx2(Te)) N*(qy1(Te)-qy2(Te))];
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

%___Plot of the normal fluxes on the interface NODES
figure(14), clf, title('Plot of the normal fluxes on interface NODES')
%Numerical plot
hold on
for i=1:length(Elements.Int)
    ielem=Elements.Int(i);
    Te=(ielem-1)*nOfElementNodes + (1:nOfElementNodes);
    N=interfaceBasisAndQuadrature.nodesdNdx(:,:,i);
    Flux1=[N*qx1(Te) N*qy1(Te)];
    Flux2=[N*qx2(Te) N*qy2(Te)];
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

%___Plot of the jump of u_star on the interface (x,jump(u))
figure(15), clf, title('Plot of the jump of u_star at interface nodes')
%Numerical plot
hold on
for i=1:length(Elements.Int)
    ielem=Elements.Int(i);
    Te=(ielem-1)*nOfNodes_star + (1:nOfNodes_star);
    N=interfaceBasisAndQuadrature_star.nodesN(:,:,i);
    u1I=N*u1_star(Te); u2I=N*u2_star(Te);
    jumpUe=u1I-u2I;
    xs=interfaceBasisAndQuadrature_star.nodesXY(:,1,i);
    ys=interfaceBasisAndQuadrature_star.nodesXY(:,2,i);
    theta=atan2(ys,xs); 
    %parche para +-pi
    k=find(abs(abs(theta)-pi)<1.e-10);if ~isempty(k)&any(abs(theta-theta(k))>pi/2), theta(k)=-theta(k); end
    plot(theta,jumpUe,'o-');
end
hold off
xlabel('theta'), ylabel('u1*-u2*')
