function plotDiscontinuosSolutionSuperParametric(fignum,X,T,u,referenceElement,Elements,LS)

%function PlotDiscontSol(fignum,u,X,T,LS,referenceElement,Elements)
% % % Check input
% % if nargin == 4
% %     nDegRef = 40;
% % end
% % 
% % % Plotting element (equal spaced points)
% % nodes = [];
% % h = 1/nDegRef;
% % for j = 0:nDegRef
% %     i = (0:nDegRef-j)'; aux = j*ones(size(i)); nodes = [nodes; [i aux]*h];
% % end
% % nodes = 2*nodes - 1; npoints = size(nodes,1);
% % 
% % % Delaunay triangulation of plotting element
% % elemTriRef = delaunayn(nodes);
% % %elemTriRef = delaunayn(nodes,{'QJ','QbB','Qc','Qx','Pp'});
% % 
% % nOfElementNodes = size(referenceElement.NodesCoord,1);
% % nOfElements = size(T,1);
% % 
% % % Compute shape functions at interpolation points
% % [N,kk1,kk2]=evaluateNodalBasisTri(nodes,referenceElement.NodesCoord,referenceElement.degree);
% % [NGeo,kk1,kk2]=evaluateNodalBasisTri(nodes,referenceElement.NodesCoordGeo,referenceElement.degreeGeo);

% Loop in elements
%patchHandle = zeros(1,nOfElements);
% % vector=[Elements.D1;Elements.D2];
% % for i = 1:length(vector)
% %     ielem=vector(i);
% %     % Interpolate solution and position at interpolation points
% %     Te = T(ielem,:);
% %     Xplot = NGeo*X(Te,:);
% %     ind = (ielem-1)*(2*nOfElementNodes)+(1:(nOfElementNodes));
% %     uplot = N*u(ind);
% % 
% %     % Plot interpolated solution in the element
% %     hold on
% %     trisurf(elemTriRef,Xplot(:,1),Xplot(:,2),uplot);
% %     hold off
% % end
% % axis equal
% % 
% % shading interp
% % 
% % vector=[Elements.Int];
% % for i = 1:length(vector)
% %     ielem=vector(i);
% %     % Interpolate solution and position at interpolation points
% %     Te = T(ielem,:);
% %     Xplot = NGeo*X(Te,:);
% %     ind = (ielem-1)*(2*nOfElementNodes)+(1:(2*nOfElementNodes));
% %     uplot = N*u(ind(1:nOfElementNodes));
% % 
% %     % Plot interpolated solution in the element
% %     hold on
% %     trisurf(elemTriRef,Xplot(:,1),Xplot(:,2),uplot);
% %     hold off
% % end
% % axis equal
% % 
% % shading interp


% % % % % % Check input
% % % % % if nargin == 4
% % % % %     nDegRef = 40;
% % % % % end
% % % % % 
% % % % % % Plotting element (equal spaced points)
% % % % % nodes = [];
% % % % % h = 1/nDegRef;
% % % % % for j = 0:nDegRef
% % % % %     i = (0:nDegRef-j)'; aux = j*ones(size(i)); nodes = [nodes; [i aux]*h];
% % % % % end
% % % % % nodes = 2*nodes-1 ; npoints = size(nodes,1);
% % % % % 
% % % % % % Delaunay triangulation of plotting element
% % % % % elemTriRef = delaunayn(nodes);
% % % % % %elemTriRef = delaunayn(nodes,{'QJ','QbB','Qc','Qx','Pp'});
% % % % % 
% % % % % nOfElementNodes = size(referenceElement.NodesCoord,1);
% % % % % nOfElements = size(T,1);
% % % % % 
% % % % % % Compute shape functions at interpolation points
% % % % % [N,kk1,kk2]=evaluateNodalBasisTri(nodes,referenceElement.NodesCoord,referenceElement.degree);
% % % % % [NGeo,kk1,kk2]=evaluateNodalBasisTri(nodes,referenceElement.NodesCoordGeo,referenceElement.degreeGeo);
% % % % % 
% % % % % % Loop in elements
% % % % % %patchHandle = zeros(1,nOfElements);
% % % % % vector=[Elements.D1;Elements.D2];
% % % % % for i = 1:length(vector)
% % % % %     ielem=vector(i);
% % % % %     % Interpolate solution and position at interpolation points
% % % % %     Te = T(ielem,:);
% % % % %     Xplot = NGeo*X(Te,:);
% % % % %     ind = (ielem-1)*(2*nOfElementNodes)+(1:(nOfElementNodes));
% % % % %     uplot = N*u(ind);
% % % % % 
% % % % %     % Plot interpolated solution in the element
% % % % %     hold on
% % % % %     trisurf(elemTriRef,Xplot(:,1),Xplot(:,2),uplot)
% % % % %     hold off
% % % % % end
% % % % % %axis equal
% % % % % 
% % % % % shading interp
% % % % % 
% % % % % %Loop in elements
% % % % % %patchHandle = zeros(1,nOfElements);
% % % % % vector1=[Elements.Int];
% % % % % for i = 1:length(vector1)
% % % % %     ielem=vector1(i);
% % % % %     % Interpolate solution and position at interpolation points
% % % % %     Te = T(ielem,:);
% % % % %     Xplot = NGeo*X(Te,:);
% % % % %     ind = (ielem-1)*(2*nOfElementNodes)+(1:(2*nOfElementNodes));
% % % % %     uplot = [N , N.*0]*u(ind);
% % % % % 
% % % % %     % Plot interpolated solution in the element
% % % % %     hold on
% % % % %     trisurf(elemTriRef,Xplot(:,1),Xplot(:,2),uplot)
% % % % %     %patchHandle(ielem) = patch('Faces',elemTriRef,'Vertices',Xplot,'FaceVertexCData',uplot,...
% % % % %     %    'FaceColor','interp','EdgeAlpha',0);
% % % % %     hold off
% % % % % end
% % % % % %axis equal
% % % % % 
% % % % % shading interp


% Input:
%  elementType: 0 for quadrilateral, 1 for triangle
elementType = 1;


p = referenceElement.degree; 
elem = elementType;
nOfElementNodes=size(referenceElement.NodesCoord,1);

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

figure(fignum); clf; hold on;
%--------------------------------------------------------------------------
% STANDARD ELEMENTS
shapeFunctions1 = computeShapeFunctionsAtPoints(p,referenceElement.NodesCoord,zpg);
N=shapeFunctions1(:,:,1)';
shapeFunctions2 =computeShapeFunctionsAtPoints(1,referenceElement.NodesCoord(1:nen_lin,:),zpg);
N_lin=shapeFunctions2(:,:,1)';

StdElements = [Elements.D1;Elements.D2];
numel = length(StdElements);
for i = 1:numel
    ielem = StdElements(i);

    Te = T(ielem,:);
    Xe = X(Te,:);
    ind = (ielem-1)*(2*nOfElementNodes) + (1:2*nOfElementNodes);
    ue = u(ind); 
    
    x = N_lin*Xe(1:nen_lin,1); 
    y = N_lin*Xe(1:nen_lin,2);
    z = N*ue(1:nOfElementNodes); 
    
    [x,y,z] = DeleteRepeated(x,y,z);
    tri = delaunay(x,y); 
    trisurf(tri,x,y,z,'FaceColor','Interp','EdgeColor','none');

end

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
    ind = (ielem-1)*(2*nOfElementNodes) + (1:2*nOfElementNodes);
    ue = u(ind);  
 [zpg,wpg,n1,n2,PtsInt] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
     zpg=zpg(1:n1,:);
    
    
    shapeFunctions3 = computeShapeFunctionsAtPoints(p,referenceElement.NodesCoord,zpg);
    N=shapeFunctions3(:,:,1)';
    N=N(1:n1,:);
    
    
    shapeFunctions4 = computeShapeFunctionsAtPoints(1,referenceElement.NodesCoord(1:nen_lin,:),zpg);
    N_lin=shapeFunctions4(:,:,1)';
    N_lin=N_lin(1:n1,:);
  
    uee=ue(1:10)+1*ue(11:end); 
    
    x = N_lin*Xe(1:nen_lin,1); 
    y = N_lin*Xe(1:nen_lin,2);
    z = N*uee;
    
    [x,y,z] = DeleteRepeated(x,y,z);
    tri = delaunay(x,y); 
    trisurf(tri,x,y,z,'FaceColor','Interp','EdgeColor','none'); 
      
end


for i = 1:numel
    ielem = EnrElements(i);

    Te = T(ielem,:);
    Xe = X(Te,:);
    LSe = LS(Te);
    ind = (ielem-1)*(2*nOfElementNodes) + (1:2*nOfElementNodes);
    ue = u(ind);  
    
     [zpg,wpg,n1,n2,PtsInt] = ModifyQuadrature(LSe,referenceElement,zgp_tri,wgp_tri,zgp_qua,wgp_qua);
     zpg=zpg(n1+1:n1+n2,:);
    
    
    shapeFunctions3 = computeShapeFunctionsAtPoints(p,referenceElement.NodesCoord,zpg);
    N=shapeFunctions3(:,:,1)';
    N=N(1:n1,:);
    
    
    shapeFunctions4 = computeShapeFunctionsAtPoints(1,referenceElement.NodesCoord(1:nen_lin,:),zpg);
    N_lin=shapeFunctions4(:,:,1)';
    N_lin=N_lin(1:n1,:);
  
    uee=ue(1:10)+-1*ue(11:end);
    x = N_lin*Xe(1:nen_lin,1); 
    y = N_lin*Xe(1:nen_lin,2);
    z = N*uee;
    
    [x,y,z] = DeleteRepeated(x,y,z);
    tri = delaunay(x,y); 
    trisurf(tri,x,y,z,'FaceColor','Interp','EdgeColor','none'); 
    

end



view(3); 
grid on; axis tight
set(gca,'FontSize',12)


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
















