function u=computeProjectionFacesSuperparametric(dataFunction,Faces,X,T,referenceElement)
%u=computeProjectionFacesSuperparametric(dataFunction,Faces,X,T,referenceEl
%ement)
% If the field contains several components, they are assummed to be given
% one variable after the other

nOfFaces = size(Faces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
faceNodes = referenceElement.faceNodes;

nOfComponents = length(feval(dataFunction,[1,1]));
nDOFFace = nOfFaceNodes*nOfComponents;

N1d = referenceElement.N1d; %Basis at approximation nodes
N1dGeo = referenceElement.N1d; Nx1dGeo = referenceElement.N1dxi; %Basis at geometry nodes
IPw_f = referenceElement.IPweights1d;
ngf = length(IPw_f);

u = zeros(nOfFaces*nOfFaceNodes*nOfComponents,1);

for f=1:nOfFaces
    iElem = Faces(f,1);  iface = Faces(f,2);
    Te = T(iElem,:);  Xe = X(Te,:);
    nodes = faceNodes(iface,:);
    xf = Xe(nodes,1);    yf = Xe(nodes,2);
    % Gauss points position
    xfg = N1dGeo*xf;  yfg = N1dGeo*yf;     
    ug = feval(dataFunction,[xfg,yfg]); 
    ug = reshape(ug,ngf,nOfComponents); %For fields with several components
    dxdxi = Nx1dGeo*xf; dydxi = Nx1dGeo*yf;
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = spdiags(dxdxiNorm.*IPw_f,0,ngf,ngf);
    M = N1d'*(dline*N1d);
    b = N1d'*(dline*ug);    
    uface = M\b;
    u(nDOFFace*(f-1)+[1:nDOFFace]) = reshape(uface,nDOFFace,1);
end
