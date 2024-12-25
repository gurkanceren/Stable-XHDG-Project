function  plotContinuosSolution(X,T,u,referenceElement)

%Check degree of the element
nOfNodesPerElement = size(T,2);

%Delaunay's solution mesh
if exist('referenceElement','var') && nOfNodesPerElement > 4
    
    % Creating solution mesh conectivity (tri)
    coordRef = referenceElement.NodesCoord;
    elemTriRef = delaunayn(coordRef);
    nOfElemTriRef = size(elemTriRef,1);
    nOfElements = size(T,1);
    tri = zeros(nOfElemTriRef*nOfElements,3);
    indexElem = 0;
    for ielem = 1:nOfElements
        Te = T(ielem,:);
        for ielemRef = 1:nOfElemTriRef
            indexElemRef = indexElem + ielemRef;
            tri(indexElemRef,:) = Te(elemTriRef(ielemRef,:));
        end
        indexElem = indexElem + nOfElemTriRef;
    end
    
else 
    tri = T;
end

%Plot
patch('Faces',tri,'Vertices',X,'FaceVertexCData',u,...
    'FaceColor','interp','EdgeAlpha',0);
axis equal
shading interp
colormap jet

