
% Academic 2D XHDG code for solving the Stokes equations with Dirichlet 
% boundary conditions. Circular interface. Dirichlet type outer bcs.
% Neumann type interface conditions. 
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
warning('off')


meshNameS = cell(4,1);
%meshNameS{1}='mesh1_P3.dcm';
meshNameS{2}='mesh1_P2.dcm';
meshNameS{3}='mesh2_P2.dcm';
meshNameS{4}='mesh3_P2.dcm';
meshNameS{5}='mesh4_P2.dcm';
errors = []; errorsPost = []; hs=[];
for i=3
    meshName=meshNameS{i};
    hs=[hs,0.5^(i+1)];
%% Load computational mesh 
%meshName = 'mesh2_P1.dcm';
if all(meshName(end-2:end)=='dcm'),GenerateMatFileFromEZ4U(['Meshes/' meshName]); end
load(['Meshes/' meshName(1:end-3) 'mat']); 

nOfElements = size(T,1); nOfElementNodesGeo = size(T,2);
X = 2*X - 1; % modify the mesh to have it defined in [-1,1]


%%Degree for approximation degreeApprox <= degreeGeo (degreeApprox==degreeGeo for isoparametric)

%% HDG preprocess
disp('XHDG preprocess...')
[F infoFaces] = hdg_preprocess(T);
nOfFaces = max(max(F)); nOfExteriorFaces = size(infoFaces.extFaces,1);
degreeApprox=2;
degreeGeo=size(elemInfo.faceNodes,2)-1;

%%reference Element
referenceElement=createReferenceElementTriSuperParametric(degreeApprox,degreeGeo);

%% Viscosity parameter: constant in every element
muElem = 1*ones(nOfElements,1);
%% Stabilization parameter
tau = repmat(1,nOfElements,3); %Uniform tau parameter (all faces)
%tau(:,2:3) = 0; %Non-null only for 1st face in each element

%% Level Sets
example = 1; % circular interface, standard element
%example = 2; %straight interface
LS = EvaluateLS(X,1);
Elements = SetElements(T,LS,[1,0],referenceElement);
figure(1),clf,
plotMesh(Elements,X,T)
figure(1), hold on;
AddFrontPlot(1,X);
hold on;
%axis tight
hold on;
x0=10;
y0=10;
width=500;
height=400;
%AddElementNumber(X,T);

 if(i==4)
%%---------------------------------------------------------------------------
% centre=[0,0];
% T1 = T(1,:); elementSize=min(norm(X(T1(1),:)-X(T1(2),:)),norm(X(T1(1),:)-X(T1(3),:)));
% [X,movedElements,movedNodes]=modifyMeshToAvoidIllConditioningCircle(centre,0.5,0.1*elementSize,X,T,referenceElement);
% LS = sqrt((X(:,1)-centre(1)).^2+(X(:,2)-centre(2)).^2)-0.5;
% %%Classification of elements for modified mesh
% Elements = SetElements(T,LS,[1,0],referenceElement);
% figure(2),clf
% plotMesh(Elements,X,T)
% figure(2); hold on;
% AddFrontPlot(example,X);
% figure(2); hold on;
% AddElementNumber(X,T);
%%------------------------------------------------------------------
 end
%% Computation
% Loop in elements
disp('Loop in elements...')

[K, f, localSolverMat, localSolverVec] = StokesHDGMatrixSuperParametric(muElem,X,T,F,referenceElement,infoFaces,tau,Elements,LS);

%Dirichlet BC
%Dirichlet face nodal coordinates
nOfFaceNodes = degreeApprox+1; nOfExteriorFaces=size(infoFaces.extFaces,1);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfExteriorFaces = size(infoFaces.extFaces,1); 
XDirichlet = zeros(nOfExteriorFaces*nOfFaceNodes,2);
if degreeApprox==degreeGeo %Isoparametric elements
    Naux=eye(nOfFaceNodes); faceNodes = referenceElement.faceNodes;
else %superparametric elements
  [Naux,kk]=evaluateNodalBasis1D(referenceElement.NodesCoord1d,referenceElement.NodesCoord1dGeo,degreeGeo); %for superparametric elements
  faceNodes = referenceElement.faceNodesGeo;
end

[uDirichlet,nullFacesint,nullFacesintindx] = computeProjectionFacesSuperparametric2(@analyticalVelocityStokes,X,T,referenceElement,Elements,LS,F,infoFaces);
dofDirichlet=[2*nOfInteriorFaces*nOfFaceNodes + [1:2*(nOfExteriorFaces)*nOfFaceNodes]];  
dofKnown = [dofDirichlet ,(2*nOfFaces*nOfFaceNodes)+[Elements.D2;Elements.Int]' , nullFacesintindx ]; uKnown =[uDirichlet ; zeros(length([Elements.D2;Elements.Int]),1);  zeros(size(nullFacesintindx,2),1)];
dofUnknown = setdiff([1:size(K,1)],dofKnown);


f = f(dofUnknown)-K(dofUnknown,dofKnown)*uKnown;
K = K(dofUnknown,dofUnknown);


% Face solution
disp('Solving linear system...')
sol = K\f;
lambda = zeros(2*nOfFaces*nOfFaceNodes+[nOfElements],1);
lambda(dofKnown)=uKnown; lambda(dofUnknown)=sol;
uhat=lambda(1:2*nOfFaces*nOfFaceNodes);
rho=[lambda(2*nOfFaces*nOfFaceNodes+1:end)];
    
% Elemental solution
disp('Calculating element by element solution...')
[u,L,p]=computeElementsSolutionStokes(uhat,rho,localSolverMat,localSolverVec,F,Elements);

% %Plot
% figure(22),clf
% plotDiscontinuosSolutionSuperParametric(X,T,u(:,1),referenceElement,20,Elements,LS)
% colorbar, title('XHDG solution: u1')
%axis tight

% figure(33),clf
% plotDiscontinuosSolutionSuperParametric(X,T,u(:,2),referenceElement,20,Elements)
% colorbar, title('XHDG solution: u2')
% axis tight

% figure(44),clf
% plotDiscontinuosSolutionSuperParametric(X,T,p,referenceElement,20,Elements)
% colorbar, title('XHDG solution: p')
% axis tight

% Local postprocess for superconvergence 
disp('Performing local postprocess...')
p_star = referenceElement.degree + 1;
nOfNodes_star = 0.5*(p_star+1)*(p_star+2);
referenceElement_star = createReferenceElementTriSuperParametric(p_star,degreeGeo);
[u_star,shapeFunctionss] = HDGpostprocessStokes(X,T,u,L,referenceElement_star,referenceElement,Elements,muElem(1),example);


% u_analytical=analyticalVelocityStokes_u1(X);
% PlotSol(5,u_analytical,X,T,LS,referenceElement,Elements)

%Plots postprocess solution
% figure(23),clf
% plotDiscontinuosSolutionSuperParametric(X,T,u_star(:,1),referenceElement_star,20,Elements)
% colorbar, title('XHDG post-solution: u1*')
% 
% figure(32),clf
% plotDiscontinuosSolutionSuperParametric(X,T,u_star(:,2),referenceElement_star,20,Elements)
% colorbar, title('XHDG post-solution: u2*')


%% Error

% analytical solution
u1_ex = @analyticalVelocityStokes_u1;
u2_ex = @analyticalVelocityStokes_u2;
p_ex=@analyticalP;

ErrorCalculations

%Error std Elements--p
% nOfElementNodes = size(referenceElement.NodesCoord,1);
% error2std = zeros(length(Elements.D1),1);
%     for i=1:length(Elements.D1)
%         iElem = [Elements.D1(i)];
%         ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
%         error2std(iElem) = computeL2NormSuperParametric(referenceElement,X(T(iElem,:),:),(1:nOfElementNodesGeo),p(ind,1),p_ex);
%     end 
% Errorstd =sqrt(sum(error2std));
% %Error Cut Elements 
% nOfElementNodes = size(referenceElement.NodesCoord,1);
% error2cut = zeros(length(Elements.Int),1);
%     for i=1:length(Elements.Int)
%         iElem = Elements.Int(i);
%         ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
%         error2cut(iElem) = computeL2NormSuperParametricCut(LS(T(iElem,:)),referenceElement,X(T(iElem,:),:),(1:nOfElementNodesGeo),p(ind,1),p_ex); 
%     end   
% Errorcut = sqrt(sum(error2cut));
% disp(['Error XHDG p= ', num2str(sqrt(Errorstd.^2+Errorcut.^2))]);
% Error_total=sqrt(Errorstd.^2+Errorcut.^2);

    
 errors = [errors, Error_total];
 errorsPost = [errorsPost, Error_total_post];
end

figure(11), clf
plot(log10(hs),log10(errors),'-o','LineWidth',1.4);
plot(log10(hs),log10(errors),'-o',log10(hs),log10(errorsPost),'o-','LineWidth',1.4);
legend('u', 'u*')
title ('P2')

slopes = (log10(errors(2:end))-log10(errors(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
slopesPost = (log10(errorsPost(2:end))-log10(errorsPost(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))








