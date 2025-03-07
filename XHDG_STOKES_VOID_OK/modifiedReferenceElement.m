

function [N,Nxi,Neta,weight,ngauss,mu,PtsInt,FaceInfo]=modifiedReferenceElement(referenceElement,LSe,mu)


p = referenceElement.degree;

%Quadrature for standart triangle and quadrilateral
[zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*p);
wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
[zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*p);

[zgp,wgp,n1,n2,PtsInt,FaceInfo] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zgp);
N = shapeFunctions(:,:,1)';  Nxi = shapeFunctions(:,:,2)'; Neta= shapeFunctions(:,:,3)';
N = N(1:n1,:);  Nxi = Nxi(1:n1,:); Neta = Neta(1:n1,:);  
weight=wgp(1:n1); ngauss=n1;    %%% ATTENTIOONN VOIID EXAMPLE 



%Enrichment of Shape Functions 


% Nl=Nold(1:n1,:)*-1;
% Nr=Nold(n1+1:n2+n1,:)*1;
% NH=[Nl;Nr];
% N=[Nold,NH];
% 
% Nlxi=Nxi_old(1:n1,:)*-1;
% Nrxi=Nxi_old(n1+1:n2+n1,:)*1;
% NHxi=[Nlxi;Nrxi];
% Nxi=[Nxi_old,NHxi];
% 
% Nleta=Neta_old(1:n1,:)*-1;
% Nreta=Neta_old(n1+1:n2+n1,:)*1;
% NHeta=[Nleta;Nreta];
% Neta=[Neta_old,NHeta];
% 
% mu=[(mu1.*ones(n1,1));(mu2*ones(n2,1))];




