function u=computeProjectionFacesSuperparametric(Faces,X,T,referenceElement,infoFaces,LS)
%u=computeProjectionFacesSuperparametric(dataFunction,Faces,X,T,referenceEl
%ement)
% If the field contains several components, they are assummed to be given
% one variable after the other

nOfFaces = size(Faces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
faceNodes = referenceElement.faceNodesGeo;
nOfInteriorFaces=size(infoFaces.intFaces,1);

N1d = referenceElement.N1d; Nx1d = referenceElement.Nxi1d; %Basis at geometry nodes
IPw_f = referenceElement.IPweights1d;
ngf = length(IPw_f);

u = zeros(nOfFaces*nOfFaceNodes*4,1);

for iface=1:nOfFaces
    elem1 = infoFaces.extFaces(iface,1); faceElem1=infoFaces.extFaces(iface,2);
    Tf = T(elem1,referenceElement.faceNodes(faceElem1,:)); Xf = X(Tf,:); LSf = LS(Tf);
   % if (all(LSf<0) | all(LSf>0)) %non-cut face
        dataFunction=@analyticalVelocityStokes;
        nOfComponents = length(feval(dataFunction,[1,1]));
        nDOFFace = nOfFaceNodes*nOfComponents;
        ind = (iface-1)*(nOfFaceNodes*4) + [1:2*nDOFFace];
        Te = T(elem1,:);  Xe = X(Te,:);
        nodes = faceNodes(faceElem1,:);
        xf = Xe(nodes,1);    yf = Xe(nodes,2);% Gauss points position
        xfg = N1d*xf;  yfg = N1d*yf;
        ug = feval(dataFunction,[xfg,yfg]);
        ug = reshape(ug,ngf,nOfComponents); %For fields with several components
        dxdxi = Nx1d*xf; dydxi = Nx1d*yf;
        dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
        dline = spdiags(dxdxiNorm.*IPw_f,0,ngf,ngf);
        M = N1d'*(dline*N1d);
        b = N1d'*(dline*ug);
        uface = M\b;
        gamaF = reshape(uface,nDOFFace,1);
        u([ind(1:nOfFaceNodes),ind(2*nOfFaceNodes+1:3*nOfFaceNodes)])=gamaF;
   
end

%Identification of zero rows (supposed to be the rows corresponding to
%enrichment in non-cut faces)
%In this implementation they are considered as Dirichlet BC with value set to 0
% nullRows = max(abs(K))<1.e-14; 
% isCCD(nullRows)=1;

%u = zeros(nOfFaces*nOfFaceNodes*nOfComponents,1);
% for f=1:nOfFaces
%     iElem = Faces(f,1);  iface = Faces(f,2);
%     Te = T(iElem,:);  Xe = X(Te,:);
%     nodes = faceNodes(iface,:);
%     xf = Xe(nodes,1);    yf = Xe(nodes,2);
%     % Gauss points position
%     xfg = N1dGeo*xf;  yfg = N1dGeo*yf;     
%     ug = feval(dataFunction,[xfg,yfg]); 
%     ug = reshape(ug,ngf,nOfComponents); %For fields with several components
%     dxdxi = Nx1dGeo*xf; dydxi = Nx1dGeo*yf;
%     dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
%     dline = spdiags(dxdxiNorm.*IPw_f,0,ngf,ngf);
%     M = N1d'*(dline*N1d);
%     b = N1d'*(dline*ug);    
%     uface = M\b;
%     u(nDOFFace*(f-1)+[1:nDOFFace]) = reshape(uface,nDOFFace,1);
% end
