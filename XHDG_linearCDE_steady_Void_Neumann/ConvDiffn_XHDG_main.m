%% NumSimsLab: 22-Mar-2023
%
%% Academic 2D X-HDG code for solving the linear Convection-Diffusion equation with 
%% Neumann boundary conditions.
%
%% Modified from 2D HDG code for Laplace eqn. with Neumann B.C. (LaCan,UPC)
%% Developed by Dr. Haroon Ahmad under the supervision of Dr. Ceren Gürkan.
%
%% This code uses the central difference scheme and the upwind scheme 
%% which could be switched using the logical flag.
%% for upwind=1 code uses upwind scheme, while for upwind=0 code uses central diff scheme.
%
% 
%
%
%
% Academic 2D HDG code for solving the Laplace equation with NEUMANN
% boundary conditions on the interface.
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

% Stabilization parameter
tau = 1;
% kinematic-viscosity:
mu = 1;
% convective-speed:
c_x = 1; 
c_y = 1;
% upwind solution or central difference
upwind = 0; %set upwind=1 for upwinding and upwind=0 for central differencing
%*****************************
errors1=[];   errorspost1=[];
%Changing mesh name (for cluster)

for p=1 %degree
    errors = []; errorsq = []; errorsPost = []; hs=[];
    for m=1:6  %mesh number
        
        filename = ['mesh' num2str(m) '_P' num2str(p) ];
        display('Solving...')
        fname1=[filename '.dcm'];
        display(fname1)
        hs=[hs,0.5^(m+1)];
        
        % Load data
        meshName = fname1;
        if all(meshName(end-2:end)=='dcm')
            GenerateMatFileFromEZ4U(['Meshes/' meshName]);
        end
        load(['Meshes/' meshName(1:end-3) 'mat']);
        %X = 2*X - 1; % modify the mesh to have it defined in [-1,1]
        
        %figure(1),clf
        %plotMesh(X,T)
        
        %Reference element
        referenceElement = createReferenceElement(1,elemInfo.nOfNodes);
        
        %X = 2*X - 1; % modify the mesh to have it defined in [-1,1]
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
        hold on;
        %AddElementNumber(X,T)
        %axis([-1,1,-1,1])
        
        
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
        [K f QQ UU Qf Uf] = hdg_matrix(X,T,F,referenceElement,infoFaces,tau,LS,Elements,mu,c_x,c_y,upwind);
        
        %System Reduction
        
        [Knew fnew uDirichlet nullFaces CD unknowns]=SystemReduction(K,f,T,X,infoFaces,referenceElement,F);
        
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
        [u_star,shapeFunctions] = hdg_postprocess(X,T,u,q,referenceElement_star,referenceElement,Elements);
        disp('Done!')

        %
        u0 = @analiticalSolutionLaplace;
        %q0 = @analiticalQLaplace;

        %Relative error for u
                
        % % Error
        %
        error_cut = zeros(length(Elements.Int),1);
        error = zeros(length(Elements.D1),1);
        
        for i = 1:length(Elements.D1) %nOfElements
            iElem=Elements.D1(i);
            ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
            error(i) = computeL2Norm(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u(ind),u0);
            
        end
        Error_D1 = sqrt(sum(error.^2));
        disp(['Error HDG D1 = ', num2str(Error_D1)]);
        
        
        for i= 1:length(Elements.Int)
            iElem=Elements.Int(i);
            ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
            error_cut(i) = computeL2Norm_cut(LS(T(iElem,:)),referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u(ind),u0);
            
        end
        Error_cut = sqrt(sum(error_cut.^2));
        disp(['Error HDG Cut = ', num2str(Error_cut)]);
        
        
        
        Error_global=sqrt(Error_cut^2+Error_D1^2);
        disp(['Error HDG for u = ', num2str(Error_global)]);
        

        %Relative error for q
                
        % % Error
        %
        %errorq_cut = zeros(length(Elements.Int),1);
        %errorq = zeros(length(Elements.D1),1);
        
        %for i = 1:length(Elements.D1) %nOfElements
        %    iElem=Elements.D1(i);
        %    ind = (iElem-1)*2*nOfElementNodes+1:iElem*2*nOfElementNodes;
        %    errorq(i) = computeL2NormQ(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),q(ind),q0);
            
        %end
        %Error_D1 = sqrt(sum(errorq.^2));
        %disp(['Error HDG D1 = ', num2str(Error_D1)]);
        
        
        %for i= 1:length(Elements.Int)
        %    iElem=Elements.Int(i);
        %    ind = (iElem-1)*2*nOfElementNodes+1:iElem*2*nOfElementNodes;
        %    errorq_cut(i) = computeL2NormQ_cut(LS(T(iElem,:)),referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),q(ind),q0);
        %    
        %end
        %Error_cut = sqrt(sum(errorq_cut.^2));
        %disp(['Error HDG Cut = ', num2str(Error_cut)]);
        
        
        
        %Errorq_global=sqrt(Error_cut^2+Error_D1^2);
        %disp(['Error HDG for q = ', num2str(Errorq_global)]);

        
        % Relative error for the postprocessed solution of u:
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
        disp(['Error HDG postprocessed D1= ', num2str(ErrorPostD1)]);

        error_postint = zeros(length(Elements.Int),1);

        for i = 1:length(Elements.Int)
            iElem=Elements.Int(i);
            Te_lin = T(iElem,:);
            Xold = X(Te_lin,:);
            Xe=shapeFunctions*Xold;
            LSe_star = EvaluateLS(Xe,1);
            ind = (iElem-1)*nOfNodes_star+1:iElem*nOfNodes_star;
            error_postint(iElem) = computeL2Norm_cut(LSe_star,referenceElement_star,Xe,(1:nOfNodes_star),u_star(ind),u0);
        end
        ErrorPostInt = sqrt(sum(error_postint.^2));
        disp(['Error HDG postprocessed Int= ', num2str(ErrorPostInt)]);
        
        Error_global_post=sqrt(ErrorPostD1^2+ErrorPostInt^2);
        disp(['Error HDG postprocessed u= ', num2str(Error_global_post)]);
        disp(' ')
        %
        % % Error
        %
        %{
        u0 = @analiticalSolutionLaplace;
        
        %Relative error in u:
        
        error_cut = zeros(length(Elements.Int),1);
        error = zeros(length(Elements.D1),1);
        
        for i = 1:length(Elements.D1) %nOfElements
            iElem=Elements.D1(i);
            ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
            error(i) = computeL2Norm(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u(ind),u0);
            
        end
        Error_D1 = sqrt(sum(error.^2));
        disp(['Error HDG D1 = ', num2str(Error_D1)]);
        
        
        for i= 1:length(Elements.Int)
            iElem=Elements.Int(i);
            ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
            error_cut(i) = computeL2Norm_cut(LS(T(iElem,:)),referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u(ind),u0);
            
        end
        Error_cut = sqrt(sum(error_cut.^2));
        disp(['Error HDG Cut = ', num2str(Error_cut)]);
        
        
        
        Error_global=sqrt(Error_cut^2+Error_D1^2);
        disp(['Error HDG  = ', num2str(Error_global)]);
        %disp('Saving the error data')
        %
        
        %Relative error in u post-processed:
        
        error_cut_post = zeros(length(Elements.Int),1);
        error_post = zeros(length(Elements.D1),1);
        
        for i = 1:length(Elements.D1) %nOfElements
            iElem=Elements.D1(i);
            ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
            error_post(i) = computeL2Norm(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u_star(ind),u0);
            
        end
        Error_D1_post = sqrt(sum(error_post.^2));
        disp(['Error HDG D1 post = ', num2str(Error_D1_post)]);
        
        
        for i= 1:length(Elements.Int)
            iElem=Elements.Int(i);
            ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
            error_cut_post(i) = computeL2Norm_cut(LS(T(iElem,:)),referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u_star(ind),u0);
            
        end
        Error_cut_post = sqrt(sum(error_cut_post.^2));
        disp(['Error HDG Cut post = ', num2str(Error_cut_post)]);
        
        
        
        Error_global_post=sqrt(Error_cut_post^2+Error_D1_post^2);
        disp(['Error HDG post  = ', num2str(Error_global_post)]);
        %disp('Saving the error data')
        %}

        %{
        %Create .output file  (for cluster)
        
        file=fopen(strcat('./', filename,'.output'),'w');
        
        fprintf(file,'Error_D1\n');
        fprintf(file, '%g\n' , Error_D1);
        fprintf(file,'Error_cut\n');
        fprintf(file, '%g\n' , Error_cut);
        fprintf(file,'Error_global\n');
        fprintf(file, '%g\n' , Error_global);
        fprintf(file,'Condition_number\n');
        fprintf(file, '%g\n' , condNu);
        
        fclose(file);
        %}
        
        
        % Plot solution
        figure,clf
        plotDiscontinuosSolution(X,T,u,referenceElement,20)
        colorbar
        title('HDG solution')
        caxis([-1 1]);
        
        figure,clf
        plotContinuosSolution(X,T,analiticalSolutionLaplace(X),referenceElement)
        colorbar
        title('HDG solution analytical')
        caxis([-1 1]);
        
        %plot faces
        %{
        vec=[infoFaces.intFaces(:,(1:2));infoFaces.extFaces];
        XDirichlet_new=computeNodalCoordinatesFaces(vec,X,T,referenceElement);
        uhat_analytical=analiticalSolutionLaplace(XDirichlet_new);
        figure
        for i=1:nOfFaces
            
            index=find(i==nullFaces);
            if ~isempty(index) i=i+1; end
            
            XDirichlet_new=computeNodalCoordinatesFaces(vec(i,:),X,T,referenceElement);
            uhat_analytical=analiticalSolutionLaplace(XDirichlet_new);
            
            plot3(XDirichlet_new(:,1),XDirichlet_new(:,2),uhat_analytical,'-go')
            grid on
            axis square
            hold on
            
            ind=(i-1)*nOfFaceNodes+(1:1:nOfFaceNodes);
            plot3(XDirichlet_new(:,1),XDirichlet_new(:,2),uhat(ind,:),'-r*')
            
        end
        %}
        %plot solution 3D
        %{
        u_analytical=analiticalSolutionLaplace(X);
        PlotSol(5,u_analytical,X,T,LS,referenceElement,Elements)
        title('HDG solution analytical')
        
        PlotDiscontSol(6,u,X,T,LS,referenceElement,Elements)
        title('HDG solution numerical')
        %}
        % convergencePlots
        
     %convergencePlots
        %% ConvergencePlots
        errors = [errors, Error_global];
        %errorsq = [errorsq, Errorq_global];
        errorsPost = [errorsPost, Error_global_post]; 
        
    end
    %% Other plots
    figure(20), clf
    plot(log10(hs),log10(errors),'-o',log10(hs),log10(errorsPost),'o-');
    legend('u','u*')

    slopes = (log10(errors(2:end))-log10(errors(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
    %slopesq = (log10(errorsq(2:end))-log10(errorsq(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
    slopesPost = (log10(errorsPost(2:end))-log10(errorsPost(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
 

    errors1=[errors1 ; errors];    
    errorspost1=[errorspost1 ; errorsPost];
    lhs = log10(hs);
end