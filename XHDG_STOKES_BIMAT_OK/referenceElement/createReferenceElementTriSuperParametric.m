function theReferenceElementSuper=createReferenceElementTriSuperParametric(nDeg,nDegGeo)

theReferenceElement=createReferenceElementTri(nDegGeo);

%Creating reference element structure
theReferenceElementSuper = struct('IPcoordinates',theReferenceElement.IPcoordinates,...
    'IPweights',theReferenceElement.IPweights,...
    'NGeo',theReferenceElement.N,...
    'NxiGeo',theReferenceElement.Nxi,...
    'NetaGeo',theReferenceElement.Neta,...
    'IPcoordinates1d',theReferenceElement.IPcoordinates1d,...
    'IPweights1d',theReferenceElement.IPweights1d',...
    'N1dGeo',theReferenceElement.N1d,'N1dxiGeo',theReferenceElement.N1dxi,...
    'faceNodesGeo',theReferenceElement.faceNodes,...
    'NodesCoordGeo',theReferenceElement.NodesCoord,...
    'NodesCoord1dGeo',theReferenceElement.NodesCoord1d,....
    'degreeGeo',nDegGeo,'degree',nDeg);

theReferenceElement=createReferenceElementTri(nDeg);

[N,Nxi,Neta]=evaluateNodalBasisTri(theReferenceElementSuper.IPcoordinates,theReferenceElement.NodesCoord,nDeg);
[N1,Nxi1]=evaluateNodalBasis1D(theReferenceElementSuper.IPcoordinates1d,theReferenceElement.NodesCoord1d,nDeg);

theReferenceElementSuper.N=N;
theReferenceElementSuper.Nxi=Nxi;
theReferenceElementSuper.Neta=Neta;
theReferenceElementSuper.N1d=N1;
theReferenceElementSuper.Nxi1d=Nxi1;
theReferenceElementSuper.NodesCoord=theReferenceElement.NodesCoord;
theReferenceElementSuper.NodesCoord1d=theReferenceElement.NodesCoord1d;
theReferenceElementSuper.faceNodes=theReferenceElement.faceNodes;