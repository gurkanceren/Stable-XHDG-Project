
isCCD = zeros(size(f)); uCCD = isCCD; 

%Identification of zero rows (supposed to be the rows corresponding to
%enrichment in non-cut faces)
%In this implementation they are considered as Dirichlet BC with value set to 0
nullRows = find(max(abs(K))<1.e-14); 
isCCD(nullRows)=1;

%Dirichlet BC
nOfInteriorFaces = size(infoFaces.intFaces,1); nOfExteriorFaces = size(infoFaces.extFaces,1);
nOfFaces=nOfInteriorFaces+nOfExteriorFaces;
nOfFaceNodes = size(referenceElement.NodesCoord1d,1); degree=referenceElement.degree;

N1d = referenceElement.N1d; %Basis at approximation nodes
IPw_f = referenceElement.IPweights1d;
Nx1d = referenceElement.N1dxi;
ngf = length(IPw_f);

for iface=1:nOfExteriorFaces
    elem1 = infoFaces.extFaces(iface,1); faceElem1=infoFaces.extFaces(iface,2);
    Tf = T(elem1,referenceElement.faceNodes(faceElem1,:)); Xf = X(Tf,:); LSf = LS(Tf);
    if (all(LSf<0) | all(LSf>0)) %non-cut face
        ind = (nOfInteriorFaces+iface-1)*(nOfFaceNodes*2) + [1:nOfFaceNodes];
        %gamaF=computeProjectionFacesSuperparametric(@analiticalSolutionLaplace,infoFaces.extFaces,X,T,referenceElement);        
        xf = Xf(:,1);    yf = Xf(:,2);
        % Gauss points position
        xfg =N1d*xf;  yfg =N1d*yf;
%         ug = feval(dataFunction,[xfg,yfg]);
%         ug = reshape(ug,ngf,nOfComponents);%For fields with several components
        ug = analiticalSolutionLaplace([xfg,yfg]);
        dxdxi = Nx1d*xf; dydxi = Nx1d*yf;
        dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
        dline = spdiags(dxdxiNorm.*IPw_f,0,ngf,ngf);
        M = N1d'*(dline*N1d);
        b = N1d'*(dline*ug);
        uface = M\b;
        gamaF=uface;
      
    else %cut face
        ind = (nOfInteriorFaces+iface-1)*(nOfFaceNodes*2) + [1:2*nOfFaceNodes];
        xf = Xf(:,1);    yf = Xf(:,2);
        % Gauss points position
        xfg =N1d*xf;  yfg =N1d*yf;
        ug = analiticalSolutionLaplace([xfg,yfg]);
        dxdxi = Nx1d*xf; dydxi = Nx1d*yf;
        dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
        dline = spdiags(dxdxiNorm.*IPw_f,0,ngf,ngf);
        M = N1d'*(dline*N1d);
        b = N1d'*(dline*ug);
        uface = M\b;
        gamaF=uface;
        gamaF = [gamaF;zeros(nOfFaceNodes,1)];
        
    end
    isCCD(ind)=1; uCCD(ind)=gamaF;
end

% System reduction (Dirichlet faces are set to prescribed value)
indCCD = find(isCCD==1); indUnknowns=setdiff([1:size(f)],indCCD);
uDirichlet = uCCD(indCCD);
f = f(indUnknowns) - K(indUnknowns,indCCD)*uDirichlet;
K=K(indUnknowns,indUnknowns);

