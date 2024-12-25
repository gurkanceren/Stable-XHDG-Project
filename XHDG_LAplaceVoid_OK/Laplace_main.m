
%% NumSimsLab: 26-Feb-2023
%% 2D XHDG code for Laplace eqn. with Dirchlet B.C. and with a Void
%% 

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

% Stabilization parameter
tau = 1;
% kinematic-viscosity:
mu = 1;

errors1=[];   errorspost1=[];
for p=1 %degree
    errors = []; errorsq = []; errorsPost = []; hs=[];
    for m=1:6 %mesh number mesh 6 and 7 are the same !! 
        
        filename = ['mesh' num2str(m) '_P' num2str(p) ];
        disp('Solving...')
        fname1=[filename '.dcm'];
        disp(fname1)
        hs=[hs,0.5^(m+1)];
        
        % Load data
        meshName = fname1;
        if all(meshName(end-2:end)=='dcm')
            GenerateMatFileFromEZ4U(['Meshes/' meshName]);
        end
        load(['Meshes/' meshName(1:end-3) 'mat']);
        X = 2*X - 1; % modify the mesh to have it defined in [-1,1]
        %figure(1),clf
        %plotMesh_old(X,T)
        %figure(1), hold on; 
        
        %Reference element
        referenceElement = createReferenceElement(1,elemInfo.nOfNodes);
        
        %Level Sets
        example = 1; % circular interface, standard element
        % example = 2; % non-standard element (l-s are two horizontal lines)
        LS = EvaluateLS(X,example);
        Elements = SetElements(T,LS,[0,1],referenceElement);
        
        figure(1),clf
        plotMesh(Elements,X,T)
        
        figure(1); hold on;
        AddFrontPlot (example);
        % % hold on;
        % %AddElementNumber(X,T)
        % axis([0,1,0,1])
         

        
        %% HDG preprocess
        disp('HDG preprocess...')
        [F, infoFaces] = hdg_preprocess(T);
        nOfElements = size(T,1);
        nOfElementNodes = size(T,2);
        nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
        nOfFaces = max(max(F));
        nOfExteriorFaces = size(infoFaces.extFaces,1);
        
        %% Computation
        % Loop in elements
        disp('Loop in elements...')
        [K, f, QQ, UU, Qf, Uf] = hdg_matrix(X,T,F,referenceElement,infoFaces,tau,LS,Elements,mu);
        
        %System Reduction
        
        [Knew, fnew, uDirichlet, CD, unknowns, nullFaces, zeroRows]=SystemReduction(K,f,T,X,infoFaces,referenceElement);
        
        condNu=condest(Knew);
        
        % Face solution
        disp('Solving linear system...')
        lambda = Knew\fnew;
        
        %Nodal Solution Rearrangement
        uhat=zeros(nOfFaceNodes*nOfFaces,1);
        uhat(unknowns)=lambda;
        uhat(CD)=uDirichlet;
        
        
        % Elemental solution
        disp('Calculating element by element solution...')
        [u,q]=computeElementsSolution(uhat,UU,QQ,Uf,Qf,T,F);

        
        %% Local postprocess for superconvergence
        % (Warning: only straight-sided elements!)
        
        disp('Performing local postprocess...')
        p_star = referenceElement.degree + 1;
        nOfNodes_star = 0.5*(p_star+1)*(p_star+2);
        referenceElement_star = createReferenceElement(1,nOfNodes_star);
        [u_star,shapeFunctions] = hdg_postprocess(X,T,u,q,referenceElement_star,referenceElement,Elements,mu);
        disp('Done!')
        
        u0 = @analiticalSolutionLaplace;

        %% Relative error for u
        error_cut = zeros(length(Elements.Int),1);
        error = zeros(length(Elements.D1),1);
        
        for i = 1:length(Elements.D1) %nOfElements
            iElem=Elements.D1(i);
            ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
            error(i) = computeL2Norm(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u(ind),u0);
            
        end
        Error_D1 = sqrt(sum(error.^2));
        %disp(['Error XHDG D1 = ', num2str(Error_D1)]);
               
        for i= 1:length(Elements.Int)
            iElem=Elements.Int(i);
            ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
            error_cut(i) = computeL2Norm_cut(LS(T(iElem,:)),referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u(ind),u0);
            
        end
        Error_cut = sqrt(sum(error_cut.^2));
        %disp(['Error XHDG Cut = ', num2str(Error_cut)]);       
        
        Error_global=sqrt(Error_cut^2+Error_D1^2);
        disp(['Error XHDG = ', num2str(Error_global)]);

        %% Relative error for the postprocessed solution of u:
        error_postd1 = zeros(length(Elements.D1),1);
        coordRef_star = referenceElement_star.NodesCoord;
        
        for i = 1:length(Elements.D1)
            iElem=Elements.D1(i);
            Te_lin = T(iElem,:);
            Xold = X(Te_lin,:);
            Xe=shapeFunctions*Xold;
            ind = (iElem-1)*nOfNodes_star+1:iElem*nOfNodes_star;
            error_postd1(iElem) = computeL2Norm(referenceElement_star,Xe,(1:nOfNodes_star),u_star(ind),u0);
        end
        ErrorPostD1 = sqrt(sum(error_postd1.^2));
        %disp(['Error HDG postprocessed D1= ', num2str(ErrorPostD1)]);

        error_postint = zeros(length(Elements.Int),1);

        for i = 1:length(Elements.Int)
            iElem=Elements.Int(i);
            Te_lin = T(iElem,:);
            Xold = X(Te_lin,:);
            Xe=shapeFunctions*Xold;
            LSe_star = EvaluateLS(Xe,example);
            ind = (iElem-1)*nOfNodes_star+1:iElem*nOfNodes_star;
            error_postint(iElem) = computeL2Norm_cut(LSe_star,referenceElement_star,Xe,(1:nOfNodes_star),u_star(ind),u0);
        end
        ErrorPostInt = sqrt(sum(error_postint.^2));
        %disp(['Error HDG postprocessed Int= ', num2str(ErrorPostInt)]);
          
        Error_global_post=sqrt(ErrorPostD1^2+ErrorPostInt^2);
        disp(['Error XHDG postprocessed u= ', num2str(Error_global_post)]);
        disp(' ')
        
        %% ConvergencePlots
        errors = [errors, Error_global];
        errorsPost = [errorsPost, Error_global_post]; 

    end
    %% Other plots
    figure(20), clf
    plot(log10(hs),log10(errors),'-o',log10(hs),log10(errorsPost),'o-');
    legend('u','u*')
    slopes = (log10(errors(2:end))-log10(errors(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))  
    slopesPost = (log10(errorsPost(2:end))-log10(errorsPost(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
    
end



