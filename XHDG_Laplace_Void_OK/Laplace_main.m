
% Academic 2D XHDG code for solving the Laplace equation with Dirichlet
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

clc, clear all, close all, %home
restoredefaultpath, setpath
% Stabilization parameter
tau = 1;

%Changing mesh name (for cluster)
errors1=[];   errorspost1=[];

for p=4 %mesh number
    errors = []; errorsPost = []; hs=[];
    for m=5 %degree 
        
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
        X = 2*X - 1; % modify the mesh to have it defined in [-1,1]
        
        
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
        axis([-1,1,-1,1])
     

        
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
        condest(K);

        %System Reduction
        [Knew fnew uDirichlet CD unknowns nullFaces zeroRows]=SystemReduction(K,f,T,X,infoFaces,referenceElement);

        % Saved to be preconditioned
        % [row, col, v] = find(Knew);
        % writematrix(fnew,'rhs.txt')
        % writematrix(row,'matrix_rows.txt')
        % writematrix(col,'matrix_cols.txt')
        % writematrix(v,'matrix_values.txt')

        %save ( 'sparse_boundary.txt','row', 'col', 'v');
        %save ('f_boundary.txt', 'fnew');

        condNu=condest(Knew);

        % Face solution
        disp('Solving linear system...')
        lambda = Knew\fnew;


        %after preconditioning
        %preconditioned = readmatrix('CG-Solution.txt');
        %preconditioned = readmatrix('CR-SOL.txt');
        %preconditioned = readmatrix('BiCGSTAB-SOL.txt');
        %preconditioned = readmatrix('ILU0-BiCGSTAB-Solution.txt');
        %preconditioned = readmatrix('ILU0-CGS-Solution.txt');
        %preconditioned = readmatrix('ILU0-CG-SOL-1D-16.txt');
        %lambda=preconditioned;

        
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
        
        u0 = @analiticalSolutionLaplace;
        
        %Relative error
                
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
        %disp(['Error HDG D1 = ', num2str(Error_D1)]);
        
        
        for i= 1:length(Elements.Int)
            iElem=Elements.Int(i);
            ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
            error_cut(i) = computeL2Norm_cut(LS(T(iElem,:)),referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u(ind),u0);
            
        end
        Error_cut = sqrt(sum(error_cut.^2));
        %disp(['Error HDG Cut = ', num2str(Error_cut)]);
        
        
        
        Error_global=sqrt(Error_cut^2+Error_D1^2);
        disp(['Error HDG  = ', num2str(Error_global)]);
        
        
        % Relative error for the postprocessed solution
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
            LSe_star = EvaluateLS(Xe,1);
            ind = (iElem-1)*nOfNodes_star+1:iElem*nOfNodes_star;
            error_postint(iElem) = computeL2Norm_cut(LSe_star,referenceElement_star,Xe,(1:nOfNodes_star),u_star(ind),u0);
        end
        ErrorPostInt = sqrt(sum(error_postint.^2));
        %disp(['Error HDG postprocessed Int= ', num2str(ErrorPostInt)]);
        
        Error_global_post=sqrt(ErrorPostD1^2+ErrorPostInt^2);
        disp(['Error HDG postprocessed = ', num2str(Error_global_post)]);
        disp(' ')
        
        
       % disp('Saving the error data')
        
%         %Create .output file  (for cluster)
%         
%         file=fopen(strcat('./', filename,'.output'),'w');
%         
%         fprintf(file,'Error_D1\n');
%         fprintf(file, '%g\n' , Error_D1);
%         fprintf(file,'Error_cut\n');
%         fprintf(file, '%g\n' , Error_cut);
%         fprintf(file,'Error_global\n');
%         fprintf(file, '%g\n' , Error_global);
%         fprintf(file,'Condition_number\n');
%         fprintf(file, '%g\n' , condNu);
%         
%         fclose(file);
        
        
        
%         % Plot solution
%         figure,clf
%         plotDiscontinuosSolution(X,T,u,referenceElement,20)
%         colorbar
%         title('HDG solution')
%         caxis([-1 1]);
%         
%         figure,clf
%         plotContinuosSolution(X,T,analiticalSolutionLaplace(X),referenceElement)
%         colorbar
%         title('HDG solution analytical')
%         caxis([-1 1]);
       
        
        %plot faces
%         vec=[infoFaces.intFaces(:,(1:2));infoFaces.extFaces];
%         XDirichlet_new=computeNodalCoordinatesFaces(vec,X,T,referenceElement);
%         uhat_analytical=analiticalSolutionLaplace(XDirichlet_new);
%         uhat_analytical(zeroRows)=0;
%         figure
%         for i=1:nOfFaces
%             
%             index=find(i==nullFaces);
%             
%             if ~isempty(index) i=i+1; end
%             
%             XDirichlet_new=computeNodalCoordinatesFaces(vec(i,:),X,T,referenceElement);
%             uhat_analytical=analiticalSolutionLaplace(XDirichlet_new);
%             
%             plot3(XDirichlet_new(:,1),XDirichlet_new(:,2),uhat_analytical,'-go')
%             grid on
%             axis square
%             hold on
%             
%             ind=(i-1)*nOfFaceNodes+(1:1:nOfFaceNodes);
%             plot3(XDirichlet_new(:,1),XDirichlet_new(:,2),uhat(ind,:),'-r*')
%             
%         end
%         
%         %plot solution 3D
%         u_analytical=analiticalSolutionLaplace(X);
%         PlotSol(5,u_analytical,X,T,LS,referenceElement,Elements)
%         title('Analytical Solution')
%         colormap(gray)
%         
%         PlotDiscontSol(6,u,X,T,LS,referenceElement,Elements)
%         title('HDG solution numerical')
        
        %convergencePlots

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
 

errors1=[errors1 ; errors ];    
errorspost1=[errorspost1 ; errorsPost];
lhs = log10(hs);

end