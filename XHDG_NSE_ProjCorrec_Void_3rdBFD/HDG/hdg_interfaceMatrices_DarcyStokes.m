function [Mn,Mt] = hdg_interfaceMatrices_DarcyStokes(infoFaces,nOfInterfaceFaces,X,T,Tcg,nOfStokesNodes,referenceElement)
%[Mn,Mt] = hdg_interfaceMatrices_DarcyStokes(Faces,X,T,Tcg,referenceElement)
% MAIN ASSUMPTION: 
%  Same type and degree for 1D approximation in the Stokes domain along 
%  the interface. Otherwise 1D shape functions for the Stokes
%  domain should be given.
%  The ordering of the nodes is given by the HDG reference element

nOfInteriorFaces = size(infoFaces.intFaces,1);
Faces = infoFaces.extFaces;
nOfFaces = nOfInteriorFaces+size(Faces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
faceNodes = referenceElement.faceNodes;
NNnx = zeros(nOfFaceNodes);
NNny = zeros(nOfFaceNodes);
NNt1t1 = zeros(nOfFaceNodes);
NNt1t2 = zeros(nOfFaceNodes);
NNt2t2 = zeros(nOfFaceNodes);
n = nOfFaceNodes*nOfInterfaceFaces;
Mn = spalloc(nOfFaceNodes*nOfFaces,2*nOfStokesNodes,2*n^2);
Mt = spalloc(2*nOfStokesNodes,2*nOfStokesNodes,4*n^2);

N1d = referenceElement.N1d; Nx1d = referenceElement.N1dxi;
IPw_f = referenceElement.IPweights1d;
ngauss_f = length(IPw_f);

for f=1:nOfInterfaceFaces
    NNnx=NNnx*0; NNny=NNny*0; NNt1t1=NNt1t1*0; NNt1t2=NNt1t2*0; NNt2t2=NNt2t2*0;
    iElem = Faces(f,1);  iface = Faces(f,2);
    Te = T(iElem,:);  Xe = X(Te,:);
    nodes = faceNodes(iface,:);
    xf = Xe(nodes,1);    yf = Xe(nodes,2);
    % Gauss points position
    xfg = N1d*xf;  yfg = N1d*yf; 
    %loop in Gauss points
    for g = 1:ngauss_f
        Nf_g = N1d(g,:);   Nfxi_g = Nx1d(g,:);
        % Integration weight
        xyDer_g = Nfxi_g*[xf yf];
        xyDerNorm_g = norm(xyDer_g);
        t = xyDer_g/xyDerNorm_g;
        nx = t(2);ny = -t(1);
        dline=IPw_f(g)*xyDerNorm_g;
        NNnx = NNnx + Nf_g'*Nf_g*(dline*nx);
        NNny = NNny + Nf_g'*Nf_g*(dline*ny);
        NNt1t1 = NNt1t1 + Nf_g'*Nf_g*(dline*t(1)^2);
        NNt2t2 = NNt2t2 + Nf_g'*Nf_g*(dline*t(2)^2);
        NNt1t2 = NNt1t2 + Nf_g'*Nf_g*(dline*t(2)*t(1));
    end
    %Assembly
    ind=nOfFaceNodes*(nOfInteriorFaces+f-1)+[1:nOfFaceNodes]; %Indexes HDG-Darcy
    Te = Tcg(f,:); ind1=2*Te-1; ind2=2*Te; %Indexes CG-velocity
    Mn(ind,ind1)=Mn(ind,ind1)+NNnx; Mn(ind,ind2)=Mn(ind,ind2)+NNny;
    Mt(ind1,ind1)=Mt(ind1,ind1)+NNt1t1;
    Mt(ind1,ind2)=Mt(ind1,ind2)+NNt1t2;
    Mt(ind2,ind2)=Mt(ind2,ind2)+NNt2t2;
    Mt(ind2,ind1)=Mt(ind1,ind2);
end