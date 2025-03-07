% LaCaN 2015 (www-lacan.upc.edu)
%
% Academic 2D X-HDG code for solving the Laplace equation BIMATERIAL PROBLEM.
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

clc, clear all, close all, %home
restoredefaultpath, setpath
warning off

% Stabilization parameter
tau = 1;
%Viscosity
global mu1 mu2 R example
mu2=100;
mu1=1;
R=0.5;

errors1=[];   errorspost1=[];

for p=3
    
%%%peraire example
meshNameS = cell(3,1);
PP=int2str(p);
meshNameS{1}=['mesh2_P',PP,'.dcm'];
meshNameS{2}=['mesh3_P',PP,'.dcm'];
meshNameS{3}=['mesh4_P',PP,'.dcm'];
meshNameS{4}=['mesh5_P',PP,'.dcm'];


errors = []; errorsPost = []; hs=[];
normsL2JumpU=[]; normsL2JumpFlux=[]; normsL2JumpUstar=[];

for i=4
    meshName=meshNameS{i};
    hs=[hs,0.5^(i+1)];
        % Load data
        if all(meshName(end-2:end)=='dcm')
            GenerateMatFileFromEZ4U(['Meshes/' meshName]);
        end
        load(['Meshes/' meshName(1:end-3) 'mat']);
        X = 2*X - 1; % modify the mesh to have it defined in [-1,1]
        
        
        %Reference element
        referenceElement = createReferenceElement(1,elemInfo.nOfNodes);
        
        %%%%%%% Start Avoid Ill conditioning %%%%%%%%%%%
        
        LS = EvaluateLS(X,1);
        Elements = SetElements(T,LS,[1,0],referenceElement);
        R=0.5;
        figure(1),clf
        plotMesh(X,T)
        hold on, theta = 0:0.01:2*pi;
        plot(R*cos(theta),R*sin(theta),'r-',X(:,1),X(:,2),'ko')
        title ('Old Mesh')
        hold off
        %---------------------------------------------------------------------------
        
%         referenceElement = createReferenceElement(1,elemInfo.nOfNodes);
%         LS = EvaluateLS(X,1);
%         Elements = SetElements(T,LS,[1,0],referenceElement);
        
        %---------------------------------------------------------------------------
        %[X,LS,movedNodes]=modifyMeshToAvoidIllConditioning(0.005,Elements,X,T,LS,referenceElement);
        [X,movedElements,movedNodes]=modifyMeshToAvoidIllConditioningCircle([0,0],R,0.1,X,T,referenceElement);
        %modifyMeshToAvoidIllConditioningCircle
        Elements = SetElements(T,LS,[1,0],referenceElement);
        figure(2),clf
        plotMesh(X,T)
        hold on, theta = 0:0.01:2*pi;
        plot(R*cos(theta),R*sin(theta),'r-',X(movedNodes,1),X(movedNodes,2),'ro')
        title ('New Mesh')
        hold off
        
        %Classification of elements for modified mesh
        Elements = SetElements(T,LS,[1,0],referenceElement);
        %%%%%%%%%%% End Avoid Ill Conditioning %%%%%%%%%%%%%%%%%%
        
        
               
        %Level Sets
        example = 1; % circular interface, standard element
        %example = 2; %straight interface
        LS = EvaluateLS(X,example);
        Elements = SetElements(T,LS,[1,0],referenceElement);
        figure(1),clf
        plotMesh(X,T)
        figure(1); hold on;
        AddFrontPlot (example,X);
        hold on;
        axis([-1 1 -1 1])
        hold on;
       % AddElementNumber(X,T);
        
        
        myvector=(1:length(X));
        non_active_nodes=setdiff(myvector,Elements.Nodes);
        
        
        %% HDG preprocess
        disp('HDG preprocess...')
        [F infoFaces] = hdg_preprocess(T);
        nOfElements = size(T,1);
        nOfElementNodes = size(T,2);
        nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
        nOfFaces = max(max(F));
        nOfExteriorFaces = size(infoFaces.extFaces,1);
        
        %% Computation
        % Loop in elements
        disp('Loop in elements...')
        [K f QQ UU Qf Uf] = x_hdg_matrix(X,T,F,referenceElement,infoFaces,tau,LS,Elements,mu1,mu2);
        
        %System Reduction
        SystemReductionScript
        
        condK=condest(K);
        
        % Face solution
        disp('Solving linear system...')
        lambda = K\f;
        
        %Nodal Solution Rearrangement
        uhat=zeros(nOfFaceNodes*2*nOfFaces,1);
        uhat(indUnknowns)=lambda;
        uhat(indCCD)=uDirichlet;
        
        % Elemental solution
        disp('Calculating element by element solution...')
        [u,q]=computeElementsSolution(uhat,UU,QQ,Uf,Qf,T,F,Elements);
       
        
        %% Local postprocess for superconvergence
        % (Warning: only straight-sided elements!)
        
        disp('Performing local postprocess...')
        p_star = referenceElement.degree + 1;
        nOfNodes_star = 0.5*(p_star+1)*(p_star+2);
        referenceElement_star = createReferenceElement(1,nOfNodes_star);
        [u_star,shapeFunctions] = hdg_postprocess4(X,T,u,q,referenceElement_star,referenceElement,Elements,mu1,mu2);
        disp('Done!')
        
        
        %% Error Calculations
        
        ErrorCalculations
        %ErrorCalculationsForQ
        
        %%Jumps on the interface: plots and integrals
        jumpsOnTheInterfaceCalculations
        normsL2JumpU=[normsL2JumpU,normL2jumpU];
        normsL2JumpFlux=[normsL2JumpFlux,normL2jumpNormalFlux]; 
        normsL2JumpUstar=[normsL2JumpUstar,normL2jumpU_star];
        plotJumpsOnInterface
        
        %% Plot solution
        PlotSolution
        
        
        %% ConvergencePlots
        
        %convergencePlots
        %convergencePlotsq
     errors = [errors, Error_total];
     errorsPost = [errorsPost, Error_total_post];       

end

figure(1), clf
plot(log10(hs),log10(errors),'-o',log10(hs),log10(errorsPost),'o-');
legend('u','u*')

slopes = (log10(errors(2:end))-log10(errors(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
slopesPost = (log10(errorsPost(2:end))-log10(errorsPost(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
 
figure(2), clf
plot(log10(hs),log10(normsL2JumpU),'-o',log10(hs),log10(normsL2JumpUstar),'o-',log10(hs),log10(normsL2JumpFlux),'o-');
legend('[u]','[u*]','[q?n]')

slopesJumpU = (log10(normsL2JumpU(2:end))-log10(normsL2JumpU(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
slopesJumpUstar = (log10(normsL2JumpUstar(2:end))-log10(normsL2JumpUstar(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
slopesJumpQn = (log10(normsL2JumpFlux(2:end))-log10(normsL2JumpFlux(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))


errors1=[errors1 ; errors ];    
errorspost1=[errorspost1 ; errorsPost];
lhs = log10(hs);

end

%converge










