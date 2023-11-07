function PlotDiscontSol(u,X,T,LS,referenceElement,Elements)

% Input:
%  elementType: 0 for quadrilateral, 1 for triangle
elementType = 1;


p = referenceElement.degree; 
elem = elementType;
nOfElementNodes=size(T,2);

if elem == 0
    nen_lin = 4;
    zVertices = referenceElement.NodesCoord;

else
    nen_lin = 3;
    zVertices = referenceElement.NodesCoord;
   
end

[zpg,wpg] = GaussLegendreCubature2D(2*p);
zpg=2*zpg-1; % mapping onto the normal reference triangle
zpg = [zpg;zVertices];


figure(1); clf; hold on; 
%--------------------------------------------------------------------------
% STANDARD ELEMENTS
shapeFunctions1 = computeShapeFunctionsAtPoints(p,referenceElement.NodesCoord,zpg);
N=shapeFunctions1(:,:,1)';
shapeFunctions2 =computeShapeFunctionsAtPoints(1,referenceElement.NodesCoord(1:nen_lin,:),zpg);
N_lin=shapeFunctions2(:,:,1)';

StdElements = [Elements.D1];
numel = length(StdElements);
for i = 1:numel
    ielem = StdElements(i);

    Te = T(ielem,:);
    Xe = X(Te,:);
    ind = (ielem-1)*nOfElementNodes+1:ielem*nOfElementNodes;
    ue = u(ind); 
    
    x = N_lin*Xe(1:nen_lin,1); 
    y = N_lin*Xe(1:nen_lin,2);
    z = N*ue; 
    
    [x,y,z] = DeleteRepeated(x,y,z);
    tri = delaunay(x,y); 
    trisurf(tri,x,y,z,'FaceColor','Interp','EdgeColor','none');

end

hold on;
%--------------------------------------------------------------------------
%CUT ELEMENTS 

[zgp_tri,wgp_tri] = GaussLegendreCubature2D(2*p);
wgp_tri = 2*wgp_tri; zgp_tri = 2*zgp_tri -1; % mapping onto the normal reference triangle
zgp_tri = [zgp_tri; 0,0; 1,0; 0,1];
wgp_tri = [wgp_tri,1,1,1];


[zgp_qua,wgp_qua] = GaussLegendreCubature2Dquad(2*p);
zgp_qua = [zgp_qua; -1,-1; 1,-1; 1,1; -1,1];
wgp_qua = [wgp_qua,1,1,1,1];


EnrElements = Elements.Int;
numel = length(EnrElements);
for i = 1:numel
    ielem = EnrElements(i);

    Te = T(ielem,:);
    Xe = X(Te,:);
    LSe = LS(Te);
    ind = (ielem-1)*nOfElementNodes+1:ielem*nOfElementNodes;
    ue = u(ind);  
    
    
    [zpg,wpg,n1,n2,PtsInt] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
     zpg=zpg(1:n1,:);
    
    
    shapeFunctions3 = computeShapeFunctionsAtPoints(p,referenceElement.NodesCoord,zpg);
    N=shapeFunctions3(:,:,1)';
    N=N(1:n1,:);
    
    
    shapeFunctions4 = computeShapeFunctionsAtPoints(1,referenceElement.NodesCoord(1:nen_lin,:),zpg);
    N_lin=shapeFunctions4(:,:,1)';
    N_lin=N_lin(1:n1,:);
  
    
    x = N_lin*Xe(1:nen_lin,1); 
    y = N_lin*Xe(1:nen_lin,2);
    z = N*ue;
    
    [x,y,z] = DeleteRepeated(x,y,z);
    tri = delaunay(x,y); 
    trisurf(tri,x,y,z,'FaceColor','Interp','EdgeColor','none');
    
    
end

view(3); 
grid on; axis tight
set(gca,'FontSize',12)
colormap("jet")


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


function [x,y,z] = DeleteRepeated(x,y,z)

npt = length(x); 
[x,ind] = sort(x); y = y(ind); z = z(ind); 
dif = abs(x(1:npt-1) - x(2:npt)) + abs(y(1:npt-1) - y(2:npt)); 
aux = find(dif < 1e-10); 
x(aux) = []; y(aux) = []; z(aux) = []; 

npt = length(y); 
[y,ind] = sort(y); x = x(ind); z = z(ind); 
dif = abs(x(1:npt-1) - x(2:npt)) + abs(y(1:npt-1) - y(2:npt)); 
aux = find(dif < 1e-10); 
x(aux) = []; y(aux) = []; z(aux) = []; 




