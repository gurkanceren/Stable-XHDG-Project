function u=computeProjectionFaces(dataFunction,Faces,X,T,referenceElement)
%should be extended to functions with several components

nOfFaces = size(Faces,1)
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
faceNodes = referenceElement.faceNodes;
M = zeros(nOfFaceNodes,nOfFaceNodes);
b = zeros(nOfFaceNodes,1);

N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;
IPw_f = referenceElement.IPweights1d;
ngauss_f = length(IPw_f);

u = [];

for f=1:nOfFaces
    M = M*0; b=b*0;
    iElem = Faces(f,1);  iface = Faces(f,2);
    Te = T(iElem,:);  Xe = X(Te,:);
    nodes = faceNodes(iface,:);
    xf = Xe(nodes,1);    yf = Xe(nodes,2);
    % Gauss points position
    xfg = N1d*xf;  yfg = N1d*yf; 
    ug = feval(dataFunction,[xfg,yfg]);
    %loop in Gauss points
    for g = 1:ngauss_f
        Nf_g = N1d(g,:);   Nfxi_g = Nx1d(g,:);
        % Integration weight
        xyDer_g = Nfxi_g*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline=IPw_f(g)*xyDerNorm_g;
        M = M + Nf_g'*Nf_g*dline;
        b = b + Nf_g'*(dline*ug(g));
    end
    
    uface = M\b;
    u = [u;uface];
end