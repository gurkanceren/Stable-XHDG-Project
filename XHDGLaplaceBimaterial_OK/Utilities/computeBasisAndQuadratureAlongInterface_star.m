function interfaceBasisAndQuadrature_star=computeBasisAndQuadratureAlongInterface_star(X,T,referenceElement,Elements,ShapeFunctions)
global example;

nDeg=referenceElement.degree;
%Quadrature for standart triangle and quadrilateral
[zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(1);
zgp_tri = referenceElement.IPcoordinates;
wgp_tri = referenceElement.IPweights;
%Quadrature for segment
N1D = referenceElement.N1d;
dNds1D = referenceElement.N1dxi;
z1D = referenceElement.IPcoordinates1d;
w1D = referenceElement.IPweights1d;
ngauss=length(z1D);

nOfCutElements = length(Elements.Int);
nOfInterfaceNodes=size(N1D,2);
nOfElementNodes=size(referenceElement.N,2);
interfaceNodesPhysical=zeros(nOfInterfaceNodes,2,nOfCutElements);
interfaceNodesRefElem=zeros(nOfInterfaceNodes,2,nOfCutElements);
interfaceIPxieta=zeros(ngauss,2,nOfCutElements);
interfaceIPxy=zeros(ngauss,2,nOfCutElements);
interfaceW=zeros(ngauss,nOfCutElements);
interfaceNnodes=zeros(nOfInterfaceNodes,nOfElementNodes,nOfCutElements);
interfacedNdxnodes=zeros(nOfInterfaceNodes,nOfElementNodes,nOfCutElements);
interfacedNdynodes=zeros(nOfInterfaceNodes,nOfElementNodes,nOfCutElements);
interfaceN=zeros(ngauss,nOfElementNodes,nOfCutElements);
interfacedNdx=zeros(ngauss,nOfElementNodes,nOfCutElements);
interfacedNdy=zeros(ngauss,nOfElementNodes,nOfCutElements);
interfaceNormal=zeros(ngauss,2,nOfCutElements);
interfaceNormalNodes=zeros(nOfInterfaceNodes,2,nOfCutElements);

%Loop in cut elements
for i=1:nOfCutElements
    Te=T(Elements.Int(i),:); %nodes in the element
    Xold=X(Te,:);
    Xe=ShapeFunctions*Xold;
    LSe=EvaluateLS(Xe,example);
    [GaussPoints,GaussWeights,n1,n2,nodesInterfaceReferenceElem,CutFaces] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
    nOfNodes = size(Xe,1);
    n=length(nodesInterfaceReferenceElem)/2; nodesInterfaceReferenceElem=reshape(nodesInterfaceReferenceElem,2,n)'; %interface nodes at reference element
    %___2D basis functions at the nodes
    [N,Nxi,Neta]=evaluateNodalBasisTri(nodesInterfaceReferenceElem,referenceElement.NodesCoord,nDeg);
    nodesInterfacePhysical= N*Xe; %interface nodes at physical element
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,n,n); invJ12 = spdiags(-J12./detJ,0,n,n);
    invJ21 = spdiags(-J21./detJ,0,n,n); invJ22 = spdiags(J11./detJ,0,n,n);
    Nx = invJ11*Nxi + invJ12*Neta; Ny = invJ21*Nxi + invJ22*Neta;
    interfaceNnodes(:,:,i)=N;
    interfacedNdxnodes(:,:,i)=Nx;
    interfacedNdynodes(:,:,i)=Ny;
    %___1D basis functions and normal vector at nodes
    [N,Nxi]=evaluateNodalBasis1D(referenceElement.NodesCoord1d,referenceElement.NodesCoord1d,nDeg);
    dxds = Nxi*nodesInterfacePhysical;
    norma = sqrt(dxds(:,1).^2+dxds(:,2).^2);
    t = dxds./[norma,norma]; n = [t(:,2),-t(:,1)];
    interfaceNormalNodes(:,:,i)=n;
    %figure(1), hold on, plot(nodesInterfacePhysical(:,1),nodesInterfacePhysical(:,2),'k*','LineWidth',2),hold off
    xieta_integrationPoints1D = N1D*nodesInterfaceReferenceElem;%coordinates of integration points in reference element
    %___2D basis functions at integration points
    [N,Nxi,Neta]=evaluateNodalBasisTri(xieta_integrationPoints1D,referenceElement.NodesCoord,nDeg);
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngauss,ngauss); invJ12 = spdiags(-J12./detJ,0,ngauss,ngauss);
    invJ21 = spdiags(-J21./detJ,0,ngauss,ngauss); invJ22 = spdiags(J11./detJ,0,ngauss,ngauss);
    Nx = invJ11*Nxi + invJ12*Neta; Ny = invJ21*Nxi + invJ22*Neta;
    xe = Xe(:,1); ye = Xe(:,2); %x and y coordinates of the element nodes
    dxds = dNds1D*nodesInterfacePhysical;
    IPxy = N1D*nodesInterfacePhysical;
    norma = sqrt(dxds(:,1).^2+dxds(:,2).^2);
    dl = w1D.*norma;
    t = dxds./[norma,norma];
    n = [t(:,2),-t(:,1)];
    
    interfaceNodesPhysical(:,:,i)=nodesInterfacePhysical;
    interfaceNodesRefElem(:,:,i)=nodesInterfaceReferenceElem;
    interfaceIPxieta(:,:,i)=xieta_integrationPoints1D;
    interfaceIPxy(:,:,i)=IPxy;
    interfaceW(:,i)=dl;
    interfaceN(:,:,i)=N;
    interfacedNdx(:,:,i)=Nx;
    interfacedNdy(:,:,i)=Ny;
    interfaceNormal(:,:,i)=n;
end

interfaceBasisAndQuadrature_star=struct('nodesXY',interfaceNodesPhysical,...
    'nodesXiEta',interfaceNodesRefElem,'IPxieta',interfaceIPxieta,...
    'IPw',interfaceW,'IPN',interfaceN,'IPdNdx',interfacedNdx,...
    'IPdNdy',interfacedNdy,'IPnormal',interfaceNormal,'IPxy',interfaceIPxy,...
    'nodesN',interfaceNnodes,'nodesdNdx',interfacedNdxnodes,...
    'nodesdNdy',interfacedNdynodes,'nodesNormal',interfaceNormalNodes);