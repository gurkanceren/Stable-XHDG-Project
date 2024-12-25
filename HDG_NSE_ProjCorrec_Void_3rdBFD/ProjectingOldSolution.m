function u=ProjectingOldSolution(dataFunction,Elements,X,T,referenceElement) %data function here is solution at old time step

nOfElements = Elements; %size(T,1);
nOfElementNodes = size(referenceElement.NodesCoord,1);
faceNodes = referenceElement.faceNodes;

nOfComponents = 1; %length(feval(dataFunction,[1,1]));
nDOFElements = nOfElementNodes*nOfComponents;

N = referenceElement.N;     
Nxi = referenceElement.Nxi; 
IPw = referenceElement.IPweights;
ng = length(IPw);

u = zeros(nOfElements*nOfElementNodes*nOfComponents,1);

for iElem=1:nOfElements
    %iElem = Elements(f,1);
    Te = T(iElem,:);  Xe = X(Te,:);
    %nodes = nOfElementNodes(iElem,:);
    xf = Xe(:,1);    yf = Xe(:,2);
    % Gauss points position
    xg = N*xf;  yg = N*yf;     
    %ug = feval(dataFunction,[xg,yg]);
    sol_old = dataFunction((iElem-1)*nOfElementNodes+(1:nOfElementNodes));
    ug = N*sol_old;
    %ug = reshape(ug,ng,nOfComponents); %For fields with several components
    dxdxi = Nxi*xf; dydxi = Nxi*yf;
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dV = spdiags(dxdxiNorm.*IPw,0,ng,ng);
    M = N'*(dV*N);
    b = N'*(dV*ug);    
    uElem = M\b;
    u(nDOFElements*(iElem-1)+[1:nDOFElements]) = reshape(uElem,nDOFElements,1);
end
