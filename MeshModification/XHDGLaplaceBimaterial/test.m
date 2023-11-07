

% A small test file to see which is the sum of the integrals over a cut element  
clc


for iElem = 3   %1:nOfElements
    
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    LSe = LS(Te);
 

c=find(iElem==Elements.D1);
d=find(iElem==Elements.D2);

if ~isempty(c)
    
    N = referenceElement.N;
    Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
    IPw = referenceElement.IPweights;
    ngauss=size(IPw,1);
    weight=IPw;
    PtsInt=[];
    FaceInfo=[];
    
    
elseif ~isempty(d)
    
    N=0;
    ngauss=0;
    PtsInt=[];
    FaceInfo=[];
    
else isempty(c) && isempty(d);
    
    
    p = referenceElement.degree;
    %Quadrature for standart triangle and quadrilateral
    [zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*p);
    wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
    [zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*p);
    
    
    [zgp,wgp,n1,n2,PtsInt,FaceInfo] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
    shapeFunctions=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord,zgp);
    N = shapeFunctions(:,:,1)';  Nxi = shapeFunctions(:,:,2)'; Neta= shapeFunctions(:,:,3)';
    ngauss1=n1; 
    ngauss2=n2;
    
end

dvolu1=0;
dvolu2=0;
% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);
    
    %% Volume computation: loop in Gauss points
for g = 1:ngauss1
      
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);    Neta_g = Neta(g,:);
    
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end
    
    %Integration weight
    dvolu1=dvolu1+wgp(g)%*det(J)
    
end
    
for g = ngauss1:(ngauss2+ngauss1)
      
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);    Neta_g = Neta(g,:);
    
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end
    
    %Integration weight
    dvolu2=dvolu2+wgp(g)  %*det(J)
  
    
end


sum=dvolu1+dvolu2



end
    
    
    
    
    
    
    