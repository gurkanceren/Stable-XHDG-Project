% Academic 2D XHDG code for solving the Stokes equations with Dirichlet 
% boundary conditions. Interface problem.
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
meshNameS{1}='mesh1_P3.dcm';
meshNameS{2}='mesh2_P3.dcm';
meshNameS{3}='mesh3_P3.dcm';
meshNameS{4}='mesh4_P3.dcm';
meshNameS{5}='mesh5_P3.dcm';
errors = []; errorsPost = []; hs=[];
for i=1:4
    meshName=meshNameS{i};
    hs=[hs,2/(2^i)];
%% Load computational mesh 
%meshName = 'mesh2_P1.dcm';
if all(meshName(end-2:end)=='dcm'),GenerateMatFileFromEZ4U(['Meshes/' meshName]); end
load(['Meshes/' meshName(1:end-3) 'mat']); 
nOfElements = size(T,1); nOfElementNodesGeo = size(T,2);
X = X*2-1 ; % modify the mesh to have it defined in [-1,1]

%%Degree for approximation <= degreeGeo (degreeApprox==degreeGeo for isoparametric)
degreeApprox=3;

%% HDG preprocess
disp('XHDG preprocess...')
[F infoFaces] = hdg_preprocess(T);
nOfFaces = max(max(F)); nOfExteriorFaces = size(infoFaces.extFaces,1);
degreeGeo=size(elemInfo.faceNodes,2)-1; 

%%reference Element
referenceElement=createReferenceElementTriSuperParametric(degreeApprox,degreeGeo);

%% Viscosity parameter
global mu1 mu2 R;
mu1 =10;
mu2 =1; 
R=0.41;
%% Stabilization parameter
tau = repmat(1,nOfElements,3); %Uniform tau parameter (all faces)
%tau(:,2:3) = 0; %Non-null only for 1st face in each element

%% Level Sets
 example = 1; % circular interface, standard element
% %example = 2; %straight interface
LS = EvaluateLS(X,example);
Elements = SetElements(T,LS,[1,0],referenceElement);
figure(1),clf,
set(gca,'FontSize',12)
x0=10;
y0=10;
width=500;
height=400;
plotMesh(X,T,Elements)
hold on;
figure(1); hold on;
AddFrontPlot (example,X);
hold on;

%u_analytical=analyticalVelocityStokes_u1(X);
% PlotSol(5,u_analytical,X,T,LS,referenceElement,Elements)

% %axis([-0.5 1.5 -0.5 1.5])
% hold on;
% %AddElementNumber(X,T);

        centre=[0,0];
        %---------------------------------------------------------------------------         
        T1 = T(1,:);
        elementSize=min(norm(X(T1(1),:)-X(T1(2),:)),norm(X(T1(1),:)-X(T1(3),:)));
        %[X,LS,movedNodes]=modifyMeshToAvoidIllConditioning(0.1*elementSize,Elements,X,T,LS,referenceElement);
        %estimate of mesh size (length of a face)
        %T1 = T(1,:); elementSize=min(norm(X(T1(1),:)-X(T1(2),:)),norm(X(T1(1),:)-X(T1(3),:)));
         [X,movedNodes]=modifyMeshToAvoidIllConditioningCircle(centre,R,0.01*elementSize,X,T,referenceElement);  
         LS = sqrt((X(:,1)-centre(1)).^2+(X(:,2)-centre(2)).^2)-R;
        
        %Classification of elements for modified mesh
        Elements = SetElements(T,LS,[1,0],referenceElement);
        
%         figure(2),clf
%         plotMesh(Elements,X,T)
%         %axis tight
%         set(gca,'FontSize',12)
%         x0=10;
%         y0=10;
%         width=500;
%         height=400;
%         set(gcf,'units','points','position',[x0,y0,width,height])
%         hold on, theta = 0:0.01:2*pi;
%         plot(R*cos(theta),R*sin(theta),'k-','LineWidth', 3)
%          hold on;
%        %AddElementNumber(X,T);
%        hold on, 
%        plot(X(:,1),X(:,2),'ko')
%        plot(X(movedNodes,1),X(movedNodes,2),'ro')
%        hold off
        %------------------------------------------------------------------



%% Computation
% Loop in elements
disp('Loop in elements...')

[K f localSolverMat localSolverVec] = StokesHDGMatrixSuperParametric(X,T,F,referenceElement,infoFaces,tau,Elements,LS);

%Dirichlet BC
%Dirichlet face nodal coordinates
nOfFaceNodes = degreeApprox+1; nOfExteriorFaces=size(infoFaces.extFaces,1);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfExteriorFaces = size(infoFaces.extFaces,1);

uDirichlet=computeProjectionFacesSuperparametric(infoFaces.extFaces,X,T,referenceElement,infoFaces,LS);
uDirichlet = [uDirichlet;0]; %mean of pressure at last element=0
dofDirichlet=[4*nOfInteriorFaces*nOfFaceNodes + [1:4*nOfExteriorFaces*nOfFaceNodes],4*nOfFaces*nOfFaceNodes+nOfElements];
dofUnknown = [[1:4*nOfInteriorFaces*nOfFaceNodes],4*nOfFaces*nOfFaceNodes+[1:nOfElements-1]];

% % System reduction (Dirichlet faces and rho at last element are set to prescribed value)
f = f(dofUnknown)-K(dofUnknown,dofDirichlet)*uDirichlet;
K = K(dofUnknown,dofUnknown);

nullRows = find(max(abs(K))<1.e-14); %enrichment in non-cut faces

% Face solution
disp('Solving linear system...')
sol = K\f;
sol(nullRows)=0; 
lambda = zeros(4*nOfFaces*nOfFaceNodes+[nOfElements],1);
lambda(dofUnknown)=sol;
lambda(dofDirichlet)=uDirichlet;
uhat=lambda(1:4*nOfFaces*nOfFaceNodes);
rho=[lambda(4*nOfFaces*nOfFaceNodes+1:end)];
    
% Elemental solution
disp('Calculating element by element solution...')
[u,L,p]=computeElementsSolutionStokes(uhat,rho,localSolverMat,localSolverVec,F,Elements);

%         %plot solution 3D
%         u_analytical=analyticalVelocityStokes_u2(X);
%         PlotSol(5,u_analytical,X,T,LS,referenceElement,Elements)
%         title('Solution Analytical')
% %         
% %         PlotDiscontSol(6,u,X,T,LS,referenceElement,Elements)
% %         title('XHDG solution')

% %Plot
% plotDiscontinuosSolutionSuperParametric(100,X,T,u(:,1),referenceElement,Elements,LS)
% colorbar, title('XHDG solution: u1')
% axis tight
% 
% plotDiscontinuosSolutionSuperParametric(101,X,T,u(:,2),referenceElement,Elements,LS)
% colorbar, title('XHDG solution: u2')
% axis tight
% 
% plotDiscontinuosSolutionSuperParametric(102,X,T,p,referenceElement,Elements,LS)
% colorbar, title('XHDG solution: p')
% axis tight

% Local postprocess for superconvergence clc
disp('Performing local postprocess...')
p_star = referenceElement.degree + 1;
nOfNodes_star = 0.5*(p_star+1)*(p_star+2);
referenceElement_star = createReferenceElementTriSuperParametric(p_star,degreeGeo);
[u_star,N] = HDGpostprocessStokes(X,T,u,L,referenceElement_star,referenceElement,Elements,example);


% %Plots postprocess solution
% plotDiscontinuosSolutionSuperParametric(103,X,T,u_star(:,1),referenceElement_star,Elements,LS)
% colorbar, title('HDG solution: u1*')
% 
% plotDiscontinuosSolutionSuperParametric(104,X,T,u_star(:,2),referenceElement_star,Elements,LS)
% colorbar, title('HDG solution: u2*')


%% Error

% analytical solution
u1_exs = @analyticalVelocityStokes_u1;
u2_exs = @analyticalVelocityStokes_u2;

%Error std Elements
nOfElementNodes = size(referenceElement.NodesCoord,1);
vector=[Elements.D1;Elements.D2];
error2std = zeros(length(vector),1);
    for i=1:length(vector)
        iElem = (vector(i));
        ind =  (iElem-1)*(2*nOfElementNodes)+(1:(2*nOfElementNodes));
        error2std(iElem) = computeL2NormSuperParametric(referenceElement,X(T(iElem,:),:),(1:nOfElementNodesGeo),u(ind,1),u1_exs).^2 ...
            +computeL2NormSuperParametric(referenceElement,X(T(iElem,:),:),(1:nOfElementNodesGeo),u(ind,2),u2_exs).^2;
    end 
Errorstd = sqrt(sum(error2std));
%Error Cut Elements 
nOfElementNodes = size(referenceElement.NodesCoord,1);
error2cut = zeros(length(Elements.Int),1);
    for i=1:length(Elements.Int)
        iElem = Elements.Int(i);
        ind =  (iElem-1)*(2*nOfElementNodes)+(1:(2*nOfElementNodes));
        error2cut(iElem) = computeL2NormSuperParametricCut(LS(T(iElem,:)),referenceElement,X(T(iElem,:),:),(1:nOfElementNodesGeo),u(ind,1),u1_exs).^2 ...
            +computeL2NormSuperParametricCut(LS(T(iElem,:)),referenceElement,X(T(iElem,:),:),(1:nOfElementNodesGeo),u(ind,2),u2_exs).^2;
    end   
Errorcut = sqrt(sum(error2cut));
disp(['Error XHDG = ', num2str(sqrt(Errorstd.^2+Errorcut.^2))]);

Error=sqrt(Errorstd.^2+Errorcut.^2);

nOfElementNodes_star = size(referenceElement_star.NodesCoord,1);
% Error for the postprocessed std Elements 
vector=[Elements.D1;Elements.D2];
error_post2 = zeros((length(vector)),1);
for i=1:length(vector);
    iElem = (vector(i));
    ind = (iElem-1)*2*nOfNodes_star+1:iElem*2*nOfNodes_star;
    Te = T(iElem,:);
    Xe = X(Te,:);
    Xe=N*Xe;
    error_post2(iElem)=computeL2NormSuperParametric(referenceElement_star,Xe,(1:nOfElementNodes_star),u_star(ind,1),u1_exs).^2 ...
    +computeL2NormSuperParametric(referenceElement_star,Xe,(1:nOfElementNodes_star),u_star(ind,2),u2_exs).^2;
end
ErrorPoststd = sqrt(sum(error_post2));

% Error for the postprocessed Cut Elements
error_post3=zeros(length(Elements.Int),1);
for i=1:length(Elements.Int)
    iElem = Elements.Int(i);
    ind = (iElem-1)*2*nOfNodes_star+1:iElem*2*nOfNodes_star;
    Te = T(iElem,:);
    Xe = X(Te,:);
    Xe=N*Xe;
    LSe_star = EvaluateLS(Xe,example);
    error_post3(iElem)=computeL2NormSuperParametricCut(LSe_star,referenceElement_star,Xe,(1:nOfElementNodes_star),u_star(ind,1),u1_exs).^2 ...
    +computeL2NormSuperParametricCut(LSe_star,referenceElement_star,Xe,(1:nOfElementNodes_star),u_star(ind,2),u2_exs).^2;
end
ErrorPostcut = sqrt(sum(error_post3));
disp(['Error XHDG postprocessed = ', num2str(sqrt(ErrorPoststd.^2+ErrorPostcut.^2))]);
ErrorPost=sqrt(ErrorPoststd.^2+ErrorPostcut.^2);
   %convergencePlots2

    errors = [errors, Error];
    errorsPost = [errorsPost, ErrorPost];
end

figure(11), clf
plot(log10(hs),log10(errors),'-o',log10(hs),log10(errorsPost),'o-','LineWidth',1.4);
legend('u','u*')


slopes = (log10(errors(2:end))-log10(errors(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
slopesPost = (log10(errorsPost(2:end))-log10(errorsPost(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))








