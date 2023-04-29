

% LaCaN 2012 (www-lacan.upc.edu)
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

%Changing mesh name (for cluster)

for m=3 %mesh number
    for p=1  %degree
        
        filename = ['mesh' num2str(m) '_P' num2str(p) ];
        display('Solving...')
        fname1=[filename '.dcm'];
        display(fname1)
        
        
        % Load data
        meshName = fname1;
        if all(meshName(end-2:end)=='dcm')
            GenerateMatFileFromEZ4U(['Meshes/' meshName]);
        end
        load(['Meshes/' meshName(1:end-3) 'mat']);
        X = 2*X - 1; % modify the mesh to have it defined in [-1,1]
        
        figure(1),clf
        plotMesh(X,T)
        
        %Reference element
        referenceElement = createReferenceElement(1,elemInfo.nOfNodes);
        
        %Level Sets
        example = 1; % circular interface, standard element
        % example = 2; % non-standard element (l-s are two horizontal lines)
        LS = EvaluateLS(X,example);
        Elements = SetElements(T,LS,[1,0],referenceElement);
        figure(1); hold on;
        AddFrontPlot (example);
        hold on;
        
        
        myvector=(1:length(X));
        non_active_nodes=setdiff(myvector,Elements.Nodes);
        
        figure(1); hold on;
        AddFrontPlot (example);
        
        
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
        [K f QQ UU Qf Uf] = hdg_matrix(X,T,F,referenceElement,infoFaces,tau,LS,Elements);
        
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
        
        % % Error
        %
        
        u0 = @analiticalSolutionLaplace;
        
        %Relative error
        
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
        disp('Saving the error data')
        
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
        
        %plot solution 3D
        u_analytical=analiticalSolutionLaplace(X);
        PlotSol(5,u_analytical,X,T,LS,referenceElement,Elements)
        title('HDG solution analytical')
        
        PlotDiscontSol(6,u,X,T,LS,referenceElement,Elements)
        title('HDG solution numerical')
        
        % convergencePlots
        
    end
end