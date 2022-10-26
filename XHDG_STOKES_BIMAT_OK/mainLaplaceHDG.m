% LaCaN 2012 (www-lacan.upc.edu) 
%
% Academic 2D HDG code for solving the Laplace equation with Dirichlet 
% boundary conditions.
%
% Main data variables:
%  X: nodal coordinates
%  T: mesh connectivitity matrix
%  F: faces (here sides) for each element 
%    (Faces are numbered so that interior faces are first)
%  elemInfo: element type
%  infoFaces.intFaces: [elem1 face1 elem2 face2 rotation] for each face
%    (Each face is assumed to have the orientation given by the 
%    first element, it is flipped when seen from the second element)
%  infoFaces.extFaces: [elem1 face1] for each exterior face
%  referenceElement: integration points and shape functions (volume and
%    boundary sides)
%

clear all, close all, %home
restoredefaultpath, setpath
    
%% Load computational mesh 
meshName = 'mesh3_P3.dcm';
if all(meshName(end-2:end)=='dcm'),GenerateMatFileFromEZ4U(['Meshes/' meshName]); end
load(['Meshes/' meshName(1:end-3) 'mat']); 
nOfElements = size(T,1); nOfElementNodesGeo = size(T,2);
figure(1),clf,plotMesh(X,T)
%% HDG preprocess
disp('HDG preprocess...')
[F infoFaces] = hdg_preprocess(T);
nOfFaces = max(max(F)); nOfExteriorFaces = size(infoFaces.extFaces,1);
degreeGeo=size(elemInfo.faceNodes,2)-1; 

%%Degree for approximation <= degreeGeo (degreeApprox==degreeGeo for isoparametric)
degreeApprox=3;

%% Diffusion parameter: constant in every element
muElem = 2*ones(nOfElements,1);
%% Stabilization parameter
tau = repmat(1,nOfElements,3); %Uniform tau parameter (all faces)
%tau(:,2:3) = 0; %Non-null only for 1st face in each element

%% Computation
% Loop in elements
disp('Loop in elements...')
if degreeApprox==degreeGeo %Isoparametric elements
    referenceElement=createReferenceElementTri(degreeApprox);
    %[K f QQ UU Qf Uf] = HDGMatrixStraightSides(muElem,X,T,F,referenceElement,infoFaces,tau);
    [K f QQ UU Qf Uf] = HDGMatrixIsoParametric(muElem,X,T,F,referenceElement,infoFaces,tau);
elseif degreeApprox<degreeGeo %Superparametric elements
    referenceElement=createReferenceElementTriSuperParametric(degreeApprox,degreeGeo);
    [K f QQ UU Qf Uf] = HDGMatrixSuperParametric(muElem,X,T,F,referenceElement,infoFaces,tau);
else
    error('degreeApprox should be <= degreeGeo')
end
%Dirichlet BC
%Dirichlet face nodal coordinates
nOfFaceNodes = degreeApprox+1; nOfExteriorFaces=size(infoFaces.extFaces,1);
nOfInteriorFaces = size(infoFaces.intFaces,1); 
XDirichlet = zeros(nOfExteriorFaces*nOfFaceNodes,2);
if degreeApprox==degreeGeo %Isoparametric elements
    Naux=eye(nOfFaceNodes); faceNodes = referenceElement.faceNodes;
else %superparametric elements
  [Naux,kk]=evaluateNodalBasis1D(referenceElement.NodesCoord1d,referenceElement.NodesCoord1dGeo,degreeGeo); %for superparametric elements
  faceNodes = referenceElement.faceNodesGeo;
end
for iFace=1:nOfExteriorFaces
    iElem = infoFaces.extFaces(iFace,1); nFace = infoFaces.extFaces(iFace,2);
    ind = (iFace-1)*nOfFaceNodes+[1:nOfFaceNodes];
    %XDirichlet(ind,:) = X(T(iElem,referenceElement.faceNodes(nFace,:)),:);
    XDirichlet(ind,:) = Naux*X(T(iElem,faceNodes(nFace,:)),:);
end
%figure(1), hold on, plot(XDirichlet(:,1),XDirichlet(:,2),'*'), hold off
uDirichlet = DirichletCondition(XDirichlet);
% System reduction (Dirichlet faces are set to prescribed value)
nDOF = nOfInteriorFaces*nOfFaceNodes;
f = f(1:nDOF)-K(1:nDOF,nDOF+1:end)*uDirichlet;
K=K(1:nDOF,1:nDOF);

% Face solution
disp('Solving linear system...')
lambda = K\f; 
uhat = [lambda; uDirichlet]; 

% Elemental solution
disp('Calculating element by element solution...')
[u,q]=computeElementsSolution(uhat,UU,QQ,Uf,Qf,F);

% Local postprocess for superconvergence 
% (Warning: only straight-sided elements!)
disp('Performing local postprocess...')
p_star = referenceElement.degree + 1;
nOfNodes_star = 0.5*(p_star+1)*(p_star+2);

if degreeApprox==degreeGeo %Isoparametric elements
    referenceElement_star = createReferenceElementTri(p_star);
    u_star = hdg_postprocess(muElem,X,T,u,q,referenceElement_star,referenceElement);
else %Superparametric elements
    referenceElement_star = createReferenceElementTriSuperParametric(p_star,degreeGeo);
    u_star = HDGpostprocessSuperParametric(muElem,X,T,u,q,referenceElement_star,referenceElement);
end
disp('Done!')

%% Error

% analytical solution
u0 = @analiticalSolutionLaplace;

nOfElementNodes = size(referenceElement.NodesCoord,1);
% Relative error
error = zeros(nOfElements,1);
if degreeApprox==degreeGeo %Isoparametric elements
    for iElem = 1:nOfElements
        ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
        error(iElem) = computeL2Norm(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u(ind),u0);
    end
else %superparametric elements
    for iElem = 1:nOfElements
        ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
        error(iElem) = computeL2NormSuperParametric(referenceElement,X(T(iElem,:),:),(1:nOfElementNodesGeo),u(ind),u0);
    end
end
    
Error = sqrt(sum(error.^2));
disp(['Error HDG = ', num2str(Error)]);

% Relative error for the postprocessed solution
error_post = zeros(nOfElements,1);

if degreeApprox==degreeGeo %Isoparametric elements
    coordRef_star = referenceElement_star.NodesCoord;
    for iElem = 1:nOfElements
        ind = (iElem-1)*nOfNodes_star+1:iElem*nOfNodes_star;
        Xe = linearMapping(X(T(iElem,1:3),:),coordRef_star); error_post(iElem) = computeL2Norm(referenceElement_star,Xe,(1:nOfNodes_star),u_star(ind),u0);
    end
else %superparametric elements
    for iElem = 1:nOfElements
        ind = (iElem-1)*nOfNodes_star+1:iElem*nOfNodes_star;
        error_post(iElem) = computeL2NormSuperParametric(referenceElement_star,X(T(iElem,:),:),(1:nOfElementNodesGeo),u_star(ind),u0);
    end
end
ErrorPost = sqrt(sum(error_post.^2));
disp(['Error HDG postprocessed = ', num2str(ErrorPost)]);
disp(' ')


%% Plot solution
if degreeApprox==degreeGeo %Isoparametric elements
    figure(2),clf
    plotDiscontinuosSolution(X,T,u,referenceElement,20)
    colorbar, title('HDG solution')
    % plot solution
    figure(3),clf
    plotPostprocessedSolution(X,T,u_star,referenceElement_star,20)
    colorbar, title('HDG postprocessed solution')
else
    figure(2),clf
    plotDiscontinuosSolutionSuperParametric(X,T,u,referenceElement,20)
    colorbar, title('HDG solution')
    % plot solution
    figure(3),clf
    plotDiscontinuosSolutionSuperParametric(X,T,u_star,referenceElement_star,20)
    colorbar, title('HDG postprocessed solution')
end




