function [u]=computeProjectionFacesSuperparametric2(dataFunction,X,T,referenceElement,Elements,LS,F,infoFaces)
%function [u,nullFaces,nullFacesintindx]=computeProjectionFacesSuperparametric2(dataFunction,X,T,referenceElement,Elements,LS,F,infoFaces)
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

for f=1:length(Faces)
    iElem = Faces(f,1);  iface = Faces(f,2);
    Te = T(iElem,:);  Xe = X(Te,:);
    nodes = faceNodes(iface,:);
    TeL=Te(nodes);
    LSe = LS(TeL);
    xf = Xe(nodes,1);    yf = Xe(nodes,2);
    
    if  LSe(1)*LSe(end)<0  &&  ismember (iElem, infoFaces.extFaces(:,1))   %% External Cut Face  
         
        [zgp_f,wgp_f,n1_f,n2_f,IntPt] = ModifyQuadrature1D(LSe,referenceElement);
        shapeFunctions_f=computeShapeFunctionsAtPoints(referenceElement.degree,referenceElement.NodesCoord1d,zgp_f);        
        ngf=n1_f;
        % Shape functions and derivatives
        N1d =shapeFunctions_f(:,:,1)';    Nx1d = shapeFunctions_f(:,:,2)';
        N1d = N1d(1:n1_f,:);     Nx1d = Nx1d(1:n1_f,:);
        %Integration weight
        IPw_f=wgp_f(1:n1_f);
        dxdxi = Nx1d*xf; dydxi = Nx1d*yf;
        
    elseif  all (xf==1)   %ismember (iElem,Elements.D2)  % external void faces
        
        continue;
        
    else  %external material face
       
        N1d = referenceElement.N1d; %Basis at approximation nodes
        Nx1d = referenceElement.Nxi1d;
        IPw_f = referenceElement.IPweights1d;
        ngf = length(IPw_f);
        dxdxi = Nx1d*xf; dydxi = Nx1d*yf;
        
 end
      
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






