


function [N,Nxi,Neta,weight,ngauss,mu,PtsInt,FaceInfo]=modifiedReferenceElement1D(referenceElement,LSe,mu1,mu2)


p = referenceElement.degree;

%Quadrature for standart triangle and quadrilateral
[zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*p);
wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
[zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*p);

[zgp,wgp,n1,n2,PtsInt,FaceInfo] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zgp);
Nold = shapeFunctions(:,:,1)';  Nxi = shapeFunctions(:,:,2)'; Neta= shapeFunctions(:,:,3)';
weight=wgp; ngauss=n1+n2;



%Enrichment of Shape Functions 

Hl=ones(n1,1);
Hr=-1*ones(n2,1);

Nl=Nold(1:n1,:)*1;
Nr=-Nold(n1+1:n2+n1,:)*1;
NH=[Nl;Nr];
N=[Nold,NH];

Nlxi=Nxi(1:n1,:)*1;
Nrxi=Nxi(n1+1:n2+n1,:)*-1;
NHxi=[Nlxi;Nrxi];
Nxi=[Nxi,NHxi];

Nleta=Neta(1:n1,:)*1;
Nreta=Neta(n1+1:n2+n1,:)*-1;
NHeta=[Nleta;Nreta];
Neta=[Neta,NHeta];

mu=[(mu1.*Hl);(mu2*Hr)];




