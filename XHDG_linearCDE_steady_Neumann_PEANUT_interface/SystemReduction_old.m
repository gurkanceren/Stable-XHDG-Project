
function [KKnew fnew uDirichlet nullFaces CD unknowns]=SystemReduction(KKold,fold,T,X,infoFaces,referenceElement,F)

%Dirichlet BC
%Dirichlet face nodal coordinates


nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfExteriorFaces = size(infoFaces.extFaces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces=size(infoFaces.intFaces,1)+size(infoFaces.extFaces,1);
degree=referenceElement.degree;

XDirichlet=computeNodalCoordinatesFaces(infoFaces.extFaces,X,T,referenceElement);
figure(1), hold on, h=plot(XDirichlet(:,1),XDirichlet(:,2),'r*');
legend([h],'Dirichlet Boundary')

%[XDirichlet,NodeDirichlet]=computeNodalCoordinatesFacesNeumann(infoFaces.extFaces,X,T,referenceElement,F);
%figure(1), hold on, plot(XDirichlet(:,1),XDirichlet(:,2),'r*')
%uDirichlet = DirichletCondition(XDirichlet);  %%%%% Dirichlet %%%%% 
%uDirichlet = NeumannCondition_2(XDirichlet);   %%%%% Neumann %%%%%% 
uDirichlet=computeProjectionFacesDirchlet(@DirichletCondition,infoFaces.extFaces,X,T,referenceElement);

% System reduction (Dirichlet faces are set to prescribed value)

%nDOF=setdiff((1:nOfFaces*nOfFaceNodes),NodeDirichlet);    %%%%% Neumann %%%%%%%
%fnew = fold(1:nDOF)-KKold(1:nDOF,nDOF+1:end)*uDirichlet;  %%%% Dirichlet
%fnew = fold(nDOF)-KKold(nDOF,NodeDirichlet)*uDirichlet;
%KKnew=KKold(nDOF,nDOF);

% System reduction (Dirichlet faces are set to prescribed value)
nDOF = nOfInteriorFaces*nOfFaceNodes;
NodeDirichlet=setdiff((1:nOfFaceNodes*nOfFaces),(1:nDOF));
fnew = fold(1:nDOF)-KKold(1:nDOF,nDOF+1:end)*uDirichlet;
KKnew=KKold(1:nDOF,1:nDOF);

%%%Further System Reduction
zeroRows=[];

%Identification of zero rows
tol=1.e-20;

Kaux=KKnew(1:degree+1:end,1:degree+1:end);
nullFaces=find(abs(sum(Kaux,1))<tol);

for i=1:length(nullFaces)
    ind = (((i)-1)*(degree+1))+(1:1:degree+1);
    zeroRows(ind)=((nullFaces(i)-1)*(degree+1))+(1:1:degree+1);
end

if ~isempty(nullFaces) 
XDirichlet_extra=computeNodalCoordinatesFaces(infoFaces.intFaces(nullFaces,1:2),X,T,referenceElement);
figure(1), hold on, plot(XDirichlet_extra(:,1),XDirichlet_extra(:,2),'c*')
end

%2nd System Reduction (Deleting zero rows and columns)

if ~isempty(zeroRows) 
KKnew(zeroRows,:)=[];
KKnew(:,zeroRows)=[];
fnew(zeroRows)=[];
end

nDOF=setdiff(nDOF,zeroRows);
CD=NodeDirichlet;
%unknowns=nDOF;
if ~isempty(zeroRows)
unknowns=setdiff((1:size(KKnew,1)),zeroRows);
else    
unknowns=(1:size(KKnew,1));    
end

disp('Hola');





























% % % function [KKnew fnew uDirichlet CD unknowns]=SystemReduction(KKold,fold,T,X,infoFaces,referenceElement)
% % % 
% % % %Dirichlet BC
% % % %Dirichlet face nodal coordinates
% % % 
% % % 
% % % nOfInteriorFaces = size(infoFaces.intFaces,1);
% % % nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
% % % nOfFaces=size(infoFaces.intFaces,1)+size(infoFaces.extFaces,1);
% % % 
% % % XDirichlet=computeNodalCoordinatesFaces(infoFaces.extFaces,X,T,referenceElement);
% % % figure(1), hold on, plot(XDirichlet(:,1),XDirichlet(:,2),'r*')
% % % uDirichlet = DirichletCondition(XDirichlet);
% % % 
% % % 
% % % % System reduction (Dirichlet faces are set to prescribed value)
% % % nDOF = nOfInteriorFaces*nOfFaceNodes;
% % % fnew = fold(1:nDOF)-KKold(1:nDOF,nDOF+1:end)*uDirichlet;
% % % KKnew=KKold(1:nDOF,1:nDOF);
% % % CD=setdiff((1:nOfFaceNodes*nOfFaces),(1:nDOF));
% % % unknowns=setdiff((1:nDOF),CD);
% % % 
% % % 
% % % % %%%Further System Reduction
% % % % zeroRows=[];
% % % % 
% % % % %Identification of zero rows
% % % % tol=1.e-20;
% % % % 
% % % % Kaux=KKnew(1:degree+1:end,1:degree+1:end);
% % % % nullFaces1=find(abs(sum(Kaux,2))<tol);
% % % % nullFaces2=find(abs(sum(Kaux,1))<tol);
% % % % nullFaces=intersect(nullFaces1,nullFaces2);
% % % % 
% % % % 
% % % % for i=1:length(nullFaces)
% % % %     ind = (((i)-1)*(degree+1))+(1:1:degree+1);
% % % %     zeroRows(ind)=((nullFaces(i)-1)*(degree+1))+(1:1:degree+1);
% % % % end
% % % % 
% % % % %Check if the identified faces are proper
% % % % 
% % % % if ~isempty(zeroRows)
% % % % unknowns=setdiff((1:size(KKnew,1)),zeroRows);
% % % % else    
% % % % unknowns=(1:size(KKnew,1));    
% % % % end
% % % % 
% % % % if ~isempty(nullFaces) 
% % % % XDirichlet_extra=computeNodalCoordinatesFaces(infoFaces.intFaces(nullFaces,1:2),X,T,referenceElement);
% % % % uDirichlet_extra = DirichletCondition(XDirichlet_extra);
% % % % figure(1), hold on, plot(XDirichlet_extra(:,1),XDirichlet_extra(:,2),'c*')
% % % % end
% % % % 
% % % % %2nd System Reduction (Deleting zero rows and columns)
% % % % 
% % % % if ~isempty(zeroRows) 
% % % % KKnew(zeroRows,:)=[];
% % % % KKnew(:,zeroRows)=[];
% % % % fnew(zeroRows)=[];
% % % % end


