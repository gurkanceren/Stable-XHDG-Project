function [u,nullFacesint,nullFacesintindx]=computeProjectionFacesSuperparametric2(dataFunction,X,T,referenceElement,Elements,LS,F,infoFaces)
%u=computeProjectionFacesSuperparametric(dataFunction,Faces,X,T,referenceEl
%ement)
% If the field contains several components, they are assummed to be given
% one variable after the other

nOfFaces = max(max(F));
Faces=infoFaces.extFaces;
FacesInt=infoFaces.intFaces;
nOfextFaces=length(Faces);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
faceNodes = referenceElement.faceNodesGeo;

nOfComponents = length(feval(dataFunction,[1,1]));
nDOFFace = nOfFaceNodes*nOfComponents;

u = zeros(nOfextFaces*nOfFaceNodes*nOfComponents,1);
nullFaces=[];
nullFacesindx =[];
nullFacesint=[];
nullFacesintindx=[];

for f=1:length(Faces)
    iElem = Faces(f,1);  iface = Faces(f,2);
    Te = T(iElem,:);  Xe = X(Te,:);
    nodes = faceNodes(iface,:);
    TeL=Te(nodes);
    LSe = LS(TeL);
    xf = Xe(nodes,1);    yf = Xe(nodes,2);
    %     figure(1); hold on;
    %     plot(xf,yf,'cd')
    
    
    N1d = referenceElement.N1d; %Basis at approximation nodes
    Nx1d = referenceElement.Nxi1d;
    IPw_f = referenceElement.IPweights1d;
    ngf = length(IPw_f);
    dxdxi = Nx1d*xf; dydxi = Nx1d*yf;
    % Gauss points position
    xfg = N1d*xf;  yfg = N1d*yf;
    ug = feval(dataFunction,[xfg,yfg]);
    ug = reshape(ug,ngf,nOfComponents); %For fields with several components
    
    
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = spdiags(dxdxiNorm.*IPw_f,0,ngf,ngf);
    M = N1d'*(dline*N1d);
    b = N1d'*(dline*ug);
    uface = M\b;
    u(nDOFFace*(f-1)+[1:nDOFFace]) = reshape(uface,nDOFFace,1);
end



%Further System Reduction 

for f=1:length(FacesInt)
    iElem1 = FacesInt(f,1);  iface1 = FacesInt(f,2);
    iElem2 = FacesInt(f,3);   iface2 = FacesInt(f,4);
    globalf1=F(iElem1,iface1);        
    Te = T(iElem1,:);  Xe = X(Te,:);
    nodes = faceNodes(iface1,:);
    xf = Xe(nodes,1);    yf = Xe(nodes,2);
    aa=ismember (iElem1,Elements.D2);
    bb=ismember (iElem2,Elements.D2);
    
    if  (aa==1 && bb==1)  % internal void faces
        nullFacesint = [nullFacesint ; [iElem1,iface1]];
        nullFacesintindx = [nullFacesintindx , nDOFFace*(globalf1-1)+[1:nDOFFace]];
%           figure (2); hold on;
%           plot(xf,yf,'rd');
              aa=[];
              bb=[];
          continue;
        
    end
              aa=[];
              bb=[];      
end


if ~isempty(nullFacesint)
    for i=1:length(nullFacesint)
Te = T(nullFacesint(i,1),:);  Xe = X(Te,:);
nodes = faceNodes(nullFacesint(i,2),:);
xf = Xe(nodes,1);    yf = Xe(nodes,2);
% figure (2); hold on;
% plot(xf,yf,'yd')
    end
end




