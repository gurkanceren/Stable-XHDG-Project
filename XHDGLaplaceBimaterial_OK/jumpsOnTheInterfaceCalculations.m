interfaceBasisAndQuadrature=computeBasisAndQuadratureAlongInterface(X,T,LS,referenceElement,Elements);
uaux=reshape(u,nOfElementNodes*2,nOfElements);
acoef=reshape(uaux(1:nOfElementNodes,:),nOfElementNodes*nOfElements,1);
bcoef=reshape(uaux(nOfElementNodes+1:end,:),nOfElementNodes*nOfElements,1);
%Standard HDG solution nodal vectors for domain 1 and 2 
u1=acoef-bcoef; u2=acoef+bcoef;

qaux=reshape(q,nOfElementNodes*4,nOfElements);
acoef_qx=reshape(qaux(1:2:2*nOfElementNodes,:),nOfElementNodes*nOfElements,1);
acoef_qy=reshape(qaux([2:2:2*nOfElementNodes],:),nOfElementNodes*nOfElements,1);
bcoef_qx=reshape(qaux(2*nOfElementNodes+[1:2:2*nOfElementNodes],:),nOfElementNodes*nOfElements,1);
bcoef_qy=reshape(qaux(2*nOfElementNodes+[2:2:2*nOfElementNodes],:),nOfElementNodes*nOfElements,1);
%Standard HDG flux vectors for domain 1 and 2 
qx1=acoef_qx-bcoef_qx; qx2=acoef_qx+bcoef_qx;
qy1=acoef_qy-bcoef_qy; qy2=acoef_qy+bcoef_qy;

%___Integral of (u1-u2)^2 and jump(normal flux)^2
integralJumpU2=0; integralJumpNormalFlux2=0;
for i=1:length(Elements.Int)
    ielem=Elements.Int(i);
    Te=(ielem-1)*nOfElementNodes + (1:nOfElementNodes);
    dl=interfaceBasisAndQuadrature.IPw(:,i);
    N=interfaceBasisAndQuadrature.IPN(:,:,i);
    u1I=N*u1(Te); u2I=N*u2(Te);
    jumpUe=u1I-u2I;
    integralJumpU2 = integralJumpU2 + dl'*(jumpUe.^2);
%    dNdx=interfaceBasisAndQuadrature.IPdNdx(:,:,i);
%    dNdy=interfaceBasisAndQuadrature.IPdNdy(:,:,i);
    jumpFlux=[N*(qx1(Te)-qx2(Te)) N*(qy1(Te)-qy2(Te))];
    nI=interfaceBasisAndQuadrature.IPnormal(:,:,i);
    jumpNormalFlux=nI(:,1).*jumpFlux(:,1)+nI(:,2).*jumpFlux(:,2);
    integralJumpNormalFlux2=integralJumpNormalFlux2+dl'*(jumpNormalFlux.^2);
end
normL2jumpU=sqrt(integralJumpU2)
normL2jumpNormalFlux=sqrt(integralJumpNormalFlux2)

%___Superconvergent solution
interfaceBasisAndQuadrature_star=computeBasisAndQuadratureAlongInterface_star(X,T,referenceElement_star,Elements,shapeFunctions);
uaux=reshape(u_star,nOfNodes_star*2,nOfElements);
acoef=reshape(uaux(1:nOfNodes_star,:),nOfNodes_star*nOfElements,1);
bcoef=reshape(uaux(nOfNodes_star+1:end,:),nOfNodes_star*nOfElements,1);
%Standard HDG solution nodal vectors for domain 1 and 2 
u1_star=acoef-bcoef; u2_star=acoef+bcoef;
%___Integral of (u1_star-u2_star)^2
integralJumpU2=0; 
for i=1:length(Elements.Int)
    ielem=Elements.Int(i);
    Te=(ielem-1)*nOfNodes_star + (1:nOfNodes_star);
    dl=interfaceBasisAndQuadrature_star.IPw(:,i);
    N=interfaceBasisAndQuadrature_star.IPN(:,:,i);
    u1I=N*u1_star(Te); u2I=N*u2_star(Te);
    jumpUe=u1I-u2I;
    integralJumpU2 = integralJumpU2 + dl'*(jumpUe.^2);
end
normL2jumpU_star=sqrt(integralJumpU2)

