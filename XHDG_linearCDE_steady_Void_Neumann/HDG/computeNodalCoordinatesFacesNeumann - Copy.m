function [XFaces,Nodes]=computeNodalCoordinatesFacesNeumann(Faces,X,T,referenceElement,F)
%XFaces=computeNodalCoordinatesFaces(Faces,X,T,referenceElement)
% 
% nOfFaces=size(Faces,1);
% nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
% XFaces = zeros(nOfFaces*nOfFaceNodes,2);
% for iFace=1:nOfFaces
%     iElem = Faces(iFace,1);
%     nFace = Faces(iFace,2);
%     ind = (iFace-1)*nOfFaceNodes+[1:nOfFaceNodes];
%     XFaces(ind,:) = X(T(iElem,referenceElement.faceNodes(nFace,:)),:);
% end


nOfFaces=size(Faces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
XFaces = zeros(nOfFaces*nOfFaceNodes,2);
Nodes= zeros(nOfFaces*nOfFaceNodes,1);
for iFace=1:nOfFaces
    
    aux=zeros(nOfFaceNodes,2);
    iElem = Faces(iFace,1);
    nFace = Faces(iFace,2);
    ind = (iFace-1)*nOfFaceNodes+[1:nOfFaceNodes];
    aux(1:nOfFaceNodes,:) = X(T(iElem,referenceElement.faceNodes(nFace,:)),:);
    
    if all((aux(1:nOfFaceNodes,1)==-1))==1 || all((aux(1:nOfFaceNodes,1)==1))==1
      
    XFaces(ind,:)=aux;
    globalfaceno=F(iElem,nFace);
    Nodes(ind,:)=(globalfaceno-1)*nOfFaceNodes+[1:nOfFaceNodes];
    
    end
        
        
end

XFaces = XFaces(any(XFaces,2),:);
Nodes(Nodes==0)=[];
%figure(1), hold on, plot(XFaces(:,1),XFaces(:,2),'*'), hold off