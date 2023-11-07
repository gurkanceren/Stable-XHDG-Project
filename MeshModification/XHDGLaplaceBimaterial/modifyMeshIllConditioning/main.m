

clc, clear all, close all, %home
restoredefaultpath, setpath


global example R;
example = 1; % circular interface, standard element
R = 0.5; 
%R = 0.567;
%R = 0.7;
%R = 0.8;
%R = 0.74;

%mesh---------------------------------------------------------------------
meshName = 'mesh5_P2.mat';
if all(meshName(end-2:end)=='dcm')
  GenerateMatFileFromEZ4U(['Meshes/' meshName]);
end
load(['Meshes/' meshName(1:end-3) 'mat']);
X = 2*X - 1; % modify the mesh to have it defined in [-1,1]

figure(1),clf
plotMesh(X,T)
hold on, theta = 0:0.01:2*pi; 
plot(R*cos(theta),R*sin(theta),'r-',X(:,1),X(:,2),'ko')
hold off
%---------------------------------------------------------------------------

referenceElement = createReferenceElement(1,elemInfo.nOfNodes);

LS = EvaluateLS(X);
Elements = SetElements(T,LS,[1,0],referenceElement); 

%---------------------------------------------------------------------------
 [X,LS,movedNodes]=modifyMeshToAvoidIllConditioning(0.1,Elements,X,T,LS,referenceElement);
 figure(2),clf
 plotMesh(X,T)
 hold on, theta = 0:0.01:2*pi; 
 plot(R*cos(theta),R*sin(theta),'r-',X(movedNodes,1),X(movedNodes,2),'ro')
 hold off

 %Classification of elements for modified mesh
 Elements = SetElements(T,LS,[1,0],referenceElement);

 

    
    







