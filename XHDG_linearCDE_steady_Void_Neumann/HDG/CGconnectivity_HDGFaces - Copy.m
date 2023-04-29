function Tcg = CGconnectivity_HDGFaces(Xcg,faces,X,T,referenceElement)
%
% Input: 
%   CG nodal coordinates: Xcg
%   HDG faces and mesh information (faces,X,T,referenceElement)
% Output:
%   Tcg: connectivity matrix for the corresponding interface/boundary
% Natural numbering is assumed for nodes in faces
nfaces = size(faces,1);
faceNodes = referenceElement.faceNodes;
nOfFaceNodes = size(faceNodes,2);
Tcg = zeros(nfaces,nOfFaceNodes);
x=Xcg(:,1); y=Xcg(:,2);
%Loop in faces
for i = 1:nfaces
    %Coordinates of nodes at i-th face
    Xf = X(T(faces(i,1),faceNodes(faces(i,2),:)),:);
    %Loop in nodes
    for j = 1:nOfFaceNodes
        [aux,k]=min( (x-Xf(j,1)).^2+(y-Xf(j,2)).^2 );
        Tcg(i,j)=k;
    end
end
%Inverse ordering for proper orientation (no if it is the second domain of an interface)
%Tcg=Tcg(:,end:-1:1); 
    