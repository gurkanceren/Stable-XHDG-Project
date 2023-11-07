function [X,movedElements,movedNodes]=modifyMeshToAvoidIllConditioningCircle(centre,radius,tol,X0,T,referenceElement)
% [Xnew,movedElements,movedNodes]=modifyMeshToAvoidIllConditioningCircle(centre,radius,tol,X,T,referenceElement)
% Modification of the mesh to avoid ill-conditioning for a circular
% interface defined by the centre and the radius
% ---> Vertex nodes closer than tol to the interface are put on it
% TO BE IMPLEMENTED
% ---> Faces are checked. If a node inside a face has different sign (including 0) than
% the 2 vertexes (face cut twice), the conditioning of the element is
% checked. If it is ill-conditioned, nodes close to the interface are
% located on top of if producing a curved face.

X = X0;
%Classsification of nodes: vertexes (1), at faces (2) and interior (0)
degree = referenceElement.degree; nOfNodes = size(X,1); nOfElements = size(T,1);
faceNodes = referenceElement.faceNodes;
interiorNodesElem = referenceElement.innerNodes;
switch degree
    case 1, nodesType = ones(nOfNodes,1);
    otherwise, nodesType = 2*ones(nOfNodes,1);
        for e=1:nOfElements
            nodesType(T(e,interiorNodesElem))=0; nodesType(T(e,[1:3]))=1;
        end
end
vertexNodes = find(nodesType==1); nodesInFaces = find(nodesType==2);

%Level-set function (signed distance to the circle)
LS = sqrt((X(:,1)-centre(1)).^2+(X(:,2)-centre(2)).^2)-radius;

%Identification of vertexes too close to the interface
nodesCloseInt = vertexNodes(find(abs(LS(vertexNodes))<tol));

%Elements cut by the interface
Elements = SetElements(T,LS,[1,0],referenceElement);
Tcut = T(Elements.Int,:);

%Loop on critical vertexes
for i=1:length(nodesCloseInt)
    iNode = nodesCloseInt(i); XiNode = X(iNode,:);
    [ii,jj]=find(Tcut==iNode); minDistance = 1.e10;
    %Loop in cut elements containing the vertex
    for k=1:length(ii);
        iElem=ii(k); iNodeElem=jj(k);
        [facesVertex,kk]=find(faceNodes==iNodeElem); %faces containing the vertex
        %Loop in faces (in cut element) containing the vertex
        for f=facesVertex'
            P1=X(Tcut(iElem,faceNodes(f,1)),:); P2=X(Tcut(iElem,faceNodes(f,end)),:);%vertexes of the face
            LS1=LS(Tcut(iElem,faceNodes(f,1))); LS2=LS(Tcut(iElem,faceNodes(f,end)));
            if(LS1*LS2<=0) %cut segment
                cutPointf = computeCutPointCircleStraightSegment(radius,centre,P1,P2);
                distancef=norm(cutPointf-XiNode);
                if(distancef<minDistance)
                    cutPoint = cutPointf; minDistance=distancef;
                end
            end
        end
    end
    if minDistance<1.e10 %it has been updated with a face :-)
        X(iNode,:)=cutPoint; LS(iNode)=0;
        %figure(1), hold on, figure(1),hold on, plot([XiNode(1),cutPoint(:,1)],[XiNode(2),cutPoint(:,2)],'r-*','LineWidth',2), hold off
    else
        error('A node has not been moved: error in the computation of the cut point in neighboring elements')
    end
end
movedNodes = nodesCloseInt;

%Elements affected by the movement
movedElements = zeros(size(T,1),1);
for k=1:length(movedNodes)
    iNode = movedNodes(k);[i,j]=find(T==iNode);movedElements(i)=1;
end
movedElements = find(movedElements==1);
movedElements1 = movedElements;

%Relocation of nodes assuming straight sides
Nlin = computeShapeFunctionsAtPoints(1,[-1,-1;1,-1;-1,1],referenceElement.NodesCoord);
Nlin = Nlin(:,:,1)';
for e=1:length(movedElements)
    Te = T(movedElements(e),:);
    Xe = X(Te(1:3),:);
    X(Te,:) = Nlin*Xe;
end

%Loop to check if there is a face with two vertexes on the interface
% -> all nodes in the face are projected to the interface
[intFaces extFaces] = GetFaces(T(:,1:3));
facesNodes = referenceElement.faceNodes;
nOfNodesFace = size(facesNodes,2);
movedElements = zeros(nOfElements,1);
for i=1:size(intFaces,1)
    ielem1=intFaces(i,1); iface1=intFaces(i,2);
    faceNodes = T(ielem1,facesNodes(iface1,:));
    nvertex1=faceNodes(1); nvertex2=faceNodes(end);%vertexes
    if LS(nvertex1)==0 & LS(nvertex2)==0
        t = X(nvertex1,:)-X(nvertex2,:); n=[t(2),-t(1)]; %not unitary!!
        for k=2:(nOfNodesFace-1)
            xk = X(faceNodes(k),:);
            cutPoint = computeCutPointCircleStraightLine(radius,centre,xk,xk+n);
            X(faceNodes(k),:)=cutPoint;
        end
        %redistribution with fekette in arclength
        X(faceNodes,:) = locatePointscurvedSegmentElement(referenceElement.NodesCoord1d,X(faceNodes,:));
        movedElements(ielem1)=1; movedElements(intFaces(i,3))=1;
        LS(faceNodes(2:end-1))=0;
    end
end
movedElements = find(movedElements==1);
movedElements = union(movedElements1,movedElements);

%Relocation of interior nodes in curved elements
XeRef = (referenceElement.NodesCoord+1)/2;
for iElem = movedElements'
    Te =T(iElem,:);
    Xe = X(Te,:);
    X(Te,:) = optimalPointsIsoparametric2D(Xe,degree,XeRef);
end

movedNodes = unique(reshape(T(movedElements,:),1,size(T,2)*length(movedElements)));
difXmovedNodes = X(movedNodes,:)-X0(movedNodes,:);
movedNodes = movedNodes(sqrt(difXmovedNodes(:,1).^2+difXmovedNodes(:,2).^2)>1.e-14);


%Loop to check other faces cut more than once
%MISSING!


%__________________________________________________________
%Function to compute the point cut by the circle in the segment P1-P2
function cutPoint = computeCutPointCircleStraightSegment(r,centre,P1,P2)

x1=P1(1); y1=P1(2);
x2=P2(1); y2=P2(2);
cx=centre(1); cy=centre(2);
alpha = (-(cx * x1) + (cx * x2) - (cy * y1) + (cy * y2) + (x1 ^ 2) - (x1 * x2) + (y1 ^ 2) - (y1 * y2) + sqrt((-cx ^ 2 * y1 ^ 2 + 2 * cx ^ 2 * y1 * y2 - cx ^ 2 * y2 ^ 2 + 2 * cx * cy * x1 * y1 - 2 * cx * cy * x1 * y2 - 2 * cx * cy * x2 * y1 + 2 * cx * cy * x2 * y2 - 2 * cx * x1 * y1 * y2 + 2 * cx * x1 * y2 ^ 2 + 2 * cx * x2 * y1 ^ 2 - 2 * cx * x2 * y1 * y2 - cy ^ 2 * x1 ^ 2 + 2 * cy ^ 2 * x1 * x2 - cy ^ 2 * x2 ^ 2 + 2 * cy * x1 ^ 2 * y2 - 2 * cy * x1 * x2 * y1 - 2 * cy * x1 * x2 * y2 + 2 * cy * x2 ^ 2 * y1 + r ^ 2 * x1 ^ 2 - 2 * r ^ 2 * x1 * x2 + r ^ 2 * x2 ^ 2 + r ^ 2 * y1 ^ 2 - 2 * r ^ 2 * y1 * y2 + r ^ 2 * y2 ^ 2 - x1 ^ 2 * y2 ^ 2 + 2 * x1 * x2 * y1 * y2 - x2 ^ 2 * y1 ^ 2))) / (x1 ^ 2 - 2 * x1 * x2 + x2 ^ 2 + y1 ^ 2 - 2 * y1 * y2 + y2 ^ 2);
cutPoint = alpha*P2+(1-alpha)*P1;
if alpha<0 | alpha>1
    alpha = -((cx * x1) - (cx * x2) + (cy * y1) - (cy * y2) - (x1 ^ 2) + (x1 * x2) - (y1 ^ 2) + (y1 * y2) + sqrt((-cx ^ 2 * y1 ^ 2 + 2 * cx ^ 2 * y1 * y2 - cx ^ 2 * y2 ^ 2 + 2 * cx * cy * x1 * y1 - 2 * cx * cy * x1 * y2 - 2 * cx * cy * x2 * y1 + 2 * cx * cy * x2 * y2 - 2 * cx * x1 * y1 * y2 + 2 * cx * x1 * y2 ^ 2 + 2 * cx * x2 * y1 ^ 2 - 2 * cx * x2 * y1 * y2 - cy ^ 2 * x1 ^ 2 + 2 * cy ^ 2 * x1 * x2 - cy ^ 2 * x2 ^ 2 + 2 * cy * x1 ^ 2 * y2 - 2 * cy * x1 * x2 * y1 - 2 * cy * x1 * x2 * y2 + 2 * cy * x2 ^ 2 * y1 + r ^ 2 * x1 ^ 2 - 2 * r ^ 2 * x1 * x2 + r ^ 2 * x2 ^ 2 + r ^ 2 * y1 ^ 2 - 2 * r ^ 2 * y1 * y2 + r ^ 2 * y2 ^ 2 - x1 ^ 2 * y2 ^ 2 + 2 * x1 * x2 * y1 * y2 - x2 ^ 2 * y1 ^ 2))) / (x1 ^ 2 - 2 * x1 * x2 + x2 ^ 2 + y1 ^ 2 - 2 * y1 * y2 + y2 ^ 2);
    cutPoint = alpha*P2+(1-alpha)*P1;
    if alpha<0 | alpha>1
        cutPoint=0;%The segment is not cut
    end
end
cutPoint = alpha*P2+(1-alpha)*P1;

%_______________________________________________________________________
%Computes the cut point in the line given by P1 and P2, closer to P1
function cutPoint = computeCutPointCircleStraightLine(r,centre,P1,P2)

x1=P1(1); y1=P1(2);
x2=P2(1); y2=P2(2);
cx=centre(1); cy=centre(2);
alpha1 = (-(cx * x1) + (cx * x2) - (cy * y1) + (cy * y2) + (x1 ^ 2) - (x1 * x2) + (y1 ^ 2) - (y1 * y2) + sqrt((-cx ^ 2 * y1 ^ 2 + 2 * cx ^ 2 * y1 * y2 - cx ^ 2 * y2 ^ 2 + 2 * cx * cy * x1 * y1 - 2 * cx * cy * x1 * y2 - 2 * cx * cy * x2 * y1 + 2 * cx * cy * x2 * y2 - 2 * cx * x1 * y1 * y2 + 2 * cx * x1 * y2 ^ 2 + 2 * cx * x2 * y1 ^ 2 - 2 * cx * x2 * y1 * y2 - cy ^ 2 * x1 ^ 2 + 2 * cy ^ 2 * x1 * x2 - cy ^ 2 * x2 ^ 2 + 2 * cy * x1 ^ 2 * y2 - 2 * cy * x1 * x2 * y1 - 2 * cy * x1 * x2 * y2 + 2 * cy * x2 ^ 2 * y1 + r ^ 2 * x1 ^ 2 - 2 * r ^ 2 * x1 * x2 + r ^ 2 * x2 ^ 2 + r ^ 2 * y1 ^ 2 - 2 * r ^ 2 * y1 * y2 + r ^ 2 * y2 ^ 2 - x1 ^ 2 * y2 ^ 2 + 2 * x1 * x2 * y1 * y2 - x2 ^ 2 * y1 ^ 2))) / (x1 ^ 2 - 2 * x1 * x2 + x2 ^ 2 + y1 ^ 2 - 2 * y1 * y2 + y2 ^ 2);
alpha2 = -((cx * x1) - (cx * x2) + (cy * y1) - (cy * y2) - (x1 ^ 2) + (x1 * x2) - (y1 ^ 2) + (y1 * y2) + sqrt((-cx ^ 2 * y1 ^ 2 + 2 * cx ^ 2 * y1 * y2 - cx ^ 2 * y2 ^ 2 + 2 * cx * cy * x1 * y1 - 2 * cx * cy * x1 * y2 - 2 * cx * cy * x2 * y1 + 2 * cx * cy * x2 * y2 - 2 * cx * x1 * y1 * y2 + 2 * cx * x1 * y2 ^ 2 + 2 * cx * x2 * y1 ^ 2 - 2 * cx * x2 * y1 * y2 - cy ^ 2 * x1 ^ 2 + 2 * cy ^ 2 * x1 * x2 - cy ^ 2 * x2 ^ 2 + 2 * cy * x1 ^ 2 * y2 - 2 * cy * x1 * x2 * y1 - 2 * cy * x1 * x2 * y2 + 2 * cy * x2 ^ 2 * y1 + r ^ 2 * x1 ^ 2 - 2 * r ^ 2 * x1 * x2 + r ^ 2 * x2 ^ 2 + r ^ 2 * y1 ^ 2 - 2 * r ^ 2 * y1 * y2 + r ^ 2 * y2 ^ 2 - x1 ^ 2 * y2 ^ 2 + 2 * x1 * x2 * y1 * y2 - x2 ^ 2 * y1 ^ 2))) / (x1 ^ 2 - 2 * x1 * x2 + x2 ^ 2 + y1 ^ 2 - 2 * y1 * y2 + y2 ^ 2);
if abs(alpha1)<abs(alpha2)
    alpha = alpha1;
else
    alpha = alpha2;
end
cutPoint = alpha*P2+(1-alpha)*P1;
