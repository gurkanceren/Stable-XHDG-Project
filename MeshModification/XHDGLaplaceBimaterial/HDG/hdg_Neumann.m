function fNeumann = hdg_Neumann(NeumannFunction,Faces,X,T,referenceElement)
%fNeumann = hdg_Neumann_laplace(@NeumannFunction,X,T,NeumannFaces)
% Output: SMALL VECTOR with size nOfNeumannFaces x nOfFaceNodes

nOfFaces = size(Faces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
faceNodes = referenceElement.faceNodes;
b = zeros(nOfFaceNodes,1);
fNeumann = zeros(nOfFaces*nOfFaceNodes,1);

N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;
IPw_f = referenceElement.IPweights1d;
ngauss_f = length(IPw_f);

for f=1:nOfFaces
    b=b*0;
    iElem = Faces(f,1);  iface = Faces(f,2);
    Te = T(iElem,:);  Xe = X(Te,:);
    nodes = faceNodes(iface,:);
    xf = Xe(nodes,1);    yf = Xe(nodes,2);
    % Gauss points position
    xfg = N1d*xf;  yfg = N1d*yf; 
    qn = feval(NeumannFunction,[xfg,yfg]);
    %loop in Gauss points
    for g = 1:ngauss_f
        Nf_g = N1d(g,:);   Nfxi_g = Nx1d(g,:);
        % Integration weight
        xyDer_g = Nfxi_g*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        dline=IPw_f(g)*xyDerNorm_g;
        b = b + Nf_g'*(dline*qn(g));
    end
    ind = (f-1)*nOfFaceNodes+[1:nOfFaceNodes];
    fNeumann(ind)=b;
end