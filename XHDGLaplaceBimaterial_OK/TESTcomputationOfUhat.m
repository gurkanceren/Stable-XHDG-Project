
allFaces = [infoFaces.intFaces(:,1:2);infoFaces.extFaces];
gama = zeros(nOfFaceNodes*2*nOfFaces,1);
for iface=1:nOfFaces
    %disp(sprintf('Face %d:',iface))
    elem1 = allFaces(iface,1); faceElem1=allFaces(iface,2);
    Tf = T(elem1,referenceElement.faceNodes(faceElem1,:)); Xf = X(Tf,:); LSf = LS(Tf);
    if (all(LSf<0) | all(LSf>0)) %non-cut face
        ind = (iface-1)*(nOfFaceNodes*2) + [1:nOfFaceNodes];
        gamaF = analiticalSolutionLaplace(Xf);
    else %cut face
        ind = (iface-1)*(nOfFaceNodes*2) + [1:2*nOfFaceNodes];
        gamaF = analyticalSolutionWithHeaviside(Xf);
    end
    gama(ind)=gamaF;
end
