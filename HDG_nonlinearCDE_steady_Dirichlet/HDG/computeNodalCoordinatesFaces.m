function XFaces=computeNodalCoordinatesFaces(Faces,X,T,referenceElement)
%XFaces=computeNodalCoordinatesFaces(Faces,X,T,referenceElement)

nOfFaces=size(Faces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
XFaces = zeros(nOfFaces*nOfFaceNodes,2);
for iFace=1:nOfFaces
    iElem = Faces(iFace,1);
    nFace = Faces(iFace,2);
    ind = (iFace-1)*nOfFaceNodes+[1:nOfFaceNodes];
    XFaces(ind,:) = X(T(iElem,referenceElement.faceNodes(nFace,:)),:);
end
%figure(1), hold on, plot(XFaces(:,1),XFaces(:,2),'*'), hold off