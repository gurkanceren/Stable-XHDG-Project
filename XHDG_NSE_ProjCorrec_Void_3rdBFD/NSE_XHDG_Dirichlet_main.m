%% NumSimLab: Sun; 13-Oct-2024
%
%% X-HDG NSE Projection Step solver with Dirichlet void in the interior.
%
%% Academic 2D X-HDG code with restart capability
%
%% Performs numerical simulation of 
%% the unsteady incompressible Navier-Stokes equations 
%% subject to Dirichlet boundary conditions for velocity at the internal void.
%% and subject to Neumann boundary conditions for correction-pressure at the internal void.
%
%% Modified from 2D X-HDG code for unsteady linear Conv-Diffn. eqn. 
%% with Dirchlet B.C. at the interior void.
%% Developed by Dr. Haroon Ahmad with Dr. Ceren Gurkan.
%% Dr. Ceren Gurkan (Asst. Prof.,Civil Engg. Dept.,Kadir Has University,Istanbul,Turkey).
%
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
%% IMPORTANT COMMENTS:
%% u is the predictor velocity 
%% uhat is the predictor velocity on the mesh skeleton
%% uc is the corrected velocity obtained in the end
%% pcorr is the correction-pressure
%% phat is the correction-pressure on the mesh skeleton
%% pc is the final pressure obtained in the end
%% Euler backward time integration
%% uses the HDG projection-Step method...
%% ... by Ueckermann & Lermusiaux (JCP,306,(2016),Pg:390-421).
%
%
clc, clear all, close all, %home
restoredefaultpath, setpath

% declare global variables
global Re
% initialize time:
time = 0;
% time-integration count:
cnt = 0;
% time-step:
dt = 1e-03;
% Reynolds number:
Re = 10;
% parameter used in discretized temporal/unsteady term:
%a_parm = 1;
% Stabilization parameter for Predictor Step:
tau = 1/Re;
% switch to write the file:
writefile = true;
% switch to restart the code after reading a previous file:
restart = false;
% maximum time:
t_max = 1.0; Time_max = 1.0;
%
%*****************************
errors1=[];   errorspost1=[];
%Changing mesh name (for cluster)

for p=1 %degree
    errors = []; errorsq = []; errorsPost = []; hs=[];
    for m=2 %mesh number 
        
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
        X(:,1) = 2*X(:,1); % modify the mesh to have it defined in [0,2]
        X(:,2) = 2*X(:,2) - 0.5;
        %figure(1),clf
        %plotMesh_old(X,T)
        %figure(1), hold on; 
        
        %Reference element
        referenceElement = createReferenceElement(1,elemInfo.nOfNodes);
        
        %Level Sets
        example = 1; % circular interface, standard element
        %example = 2; % flower/peanut interface
        LS = EvaluateLS(X,example);
        Elements = SetElements(T,LS,[0,1],referenceElement);
        %
        %% Plot Mesh 
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

        % initialize the velocity and pressure field:
        u_crt = zeros(nOfElements*nOfElementNodes,2);
        u_old = zeros(nOfElements*nOfElementNodes,2);
        u_old2 = zeros(nOfElements*nOfElementNodes,2);
        u_old3 = zeros(nOfElements*nOfElementNodes,2); 
        p_crt = zeros(nOfElements*nOfElementNodes,1);

        %% Computation
        %
        time=0;
        %
        %% Read the file for time here
        %
        if(restart==true)
            %
            %% Read the time here
            %
            fileID = fopen('timeANDcount.txt','r');
            formatSpec = '%f %d %d %d';
            sizeA = [4 Inf];
            oyge = fscanf(fileID,formatSpec,sizeA);
            fclose(fileID);
            %
            time = oyge(1);
            %
            %% Read the file for 2D flowfield here
            %
            fileID = fopen('2DFlowField.txt','r');
            formatSpec = '%e %e %e';
            sizeA = [3 Inf];
            A = fscanf(fileID,formatSpec,sizeA);
            fclose(fileID);
            %
            At = A';
            %
            u = At(:,1);
            v = At(:,2);
            p_crt = At(:,3);
            %
            u_crt(:,1) = u;
            u_crt(:,2) = v;
            %
        end
        %% Always initialize the count as zero
        cnt = 0;
         %% the time-loop starts here:
         while(time<t_max) 
         %while(cnt<100)

             %
             time = time+dt; % update time
             cnt = cnt + 1;
             %
             %
             % save previous time data
             %u_old = u_crt;
             p_old = p_crt;
             %storage of u values of old times
             if(cnt==1)      %for n=1
                 %
                 u_old = u_crt;  
                 %here: u_old = u_(n-1)
                 % multiplier for discretized unsteady velocity term at n time step
                 a_parm = 1;     
             elseif(cnt==2)  %for n=2
                 %
                 u_old2 = u_crt; 
                 %here: u_old2 = u_(n-1) and u_old = u_(n-2)
                 % multiplier for discretized unsteady velocity term at n time step
                 a_parm = 2/3;       
             elseif(cnt==3)  %for n=3
                 %
                 u_old3 = u_crt; 
                 %here: u_old3 = u_(n-1) and u_old2 = u_(n-2) and u_old = u_(n-3)
                 % multiplier for discretized unsteady velocity term at n time step
                 a_parm = 6/11;        
             else            %for n>3
                 %
                 u_old = u_old2;
                 u_old2 = u_old3;
                 u_old3 = u_crt;
                 %here: u_old3 = u_(n-1) and u_old2 = u_(n-2) and u_old = u_(n-3)
                 % multiplier for discretized unsteady velocity term at n time step
                 a_parm = 6/11;       
             end
             %
             % Stabilization parameter for Pressure Correction Step:
             tauP = 1/(a_parm*dt*tau);  
             %
             %
             %% Computation Step-1 (Predictor Step: gives Predictor velocity):
             % Loop in elementsnOfFaces
             % disp('Predictor Step: Loop in elements...')
             %
             [K f QQ UU Qf Uf] = hdg_matrix_predictorStepDirichlet(X,T,F,referenceElement,infoFaces,tau,LS,Elements,dt,cnt,u_old,u_old2,u_old3,p_old,time,Re,a_parm);
             %
             %System Reduction
             %
             [Knew fnew uDirichlet CD unknowns nullFaces zeroRows]=SystemReductionPredictorStep(K,f,T,X,infoFaces,referenceElement,time,Re);
             %
             condNu=condest(Knew);
             %
             % Face solution
             % disp('Solving linear system for predictor velocity: u...')
             lambda = Knew\fnew;
             %    
             %Nodal Solution Rearrangement
             uhat=zeros(2*nOfFaceNodes*nOfFaces,1);
             uhat(unknowns)=lambda;
             uhat(CD)=uDirichlet;
             %
             % Elemental solution
             % disp('Calculating element by element solution for predictor velocity...')
             [u,q]=computeElementsSolutionPredictorStep(uhat,UU,QQ,Uf,Qf,T,F,referenceElement);
             %               
             %% Computation Step-2 (Correction-Pressure Poisson eqn.: gives correction-pressure):
             % Loop in elementsnOfFaces
             % disp('Correction Pressure Step: Loop in elements...')
             %
             [Kp fp QQp UUp Qfp Ufp] = hdg_matrix_PressCorrNeumann(X,T,F,referenceElement,infoFaces,tauP,LS,Elements,dt,Re,a_parm,u,uhat,time);
             %
             %System Reduction
             %
             [Kpnew fpnew pDirichlet CDp unknownsp nullFacesp zeroRowsp]=SystemReductionPCorr(Kp,fp,T,X,infoFaces,referenceElement,time,Re);
             %
             condNu=condest(Kpnew);
             %
             % Face solution
             %disp('Solving linear system...')
             lambda_pcorr = Kpnew\fpnew;
             %
             %Nodal Solution Rearrangement
             phat=zeros(nOfFaceNodes*nOfFaces,1);
             phat(unknownsp)=lambda_pcorr;
             phat(CDp)=pDirichlet;
             %
             % Elemental solution
             %    disp('Calculating element by element solution...')
             [pcorr,qcorr]=computeElementsSolutionPressCorr(phat,UUp,QQp,Ufp,Qfp,T,F);
             %
             %% Computation Step-3 (Corrector Step: gives corrected velocity and corrected pressure)
             %
             % Element-by-Element correction
             %    disp('Perform element by element correction of velocity & pressure...')
             %%%[u_crt, p_crt]=computeElementalCorrectionStep(u,uhat,pcorr,phat,qcorr,X,T,F,referenceElement,infoFaces,tauP,dt,Re,a_parm,time,cnt,p_old);
             %
             %% eqn.(54) of the paper for obtaining corrected pressure gives better results
             %[u_crt,p_crt,p_wrc_crt]=computeElementalCorrectionStep_XHDG(LS,Elements,u,uhat,pcorr,phat,qcorr,X,T,F,referenceElement,infoFaces,tauP,dt,Re,a_parm,time,cnt,p_old);
             [u_crt,p_crt,p_wrc_crt]=computeElementalCorrectionStep_HDG(u,uhat,pcorr,phat,qcorr,X,T,F,referenceElement,infoFaces,tauP,dt,Re,a_parm,time,cnt,p_old);

             disp('count of time-integration')
             disp(cnt)
             %
         end  
         %% the time-loop ends here.
         %disp('Hola')
         %
         %% save the flow-field (u_corrected and p_corrected) in a text file
         %
         writefile = true;
         %
         if(writefile==true)
             %
             %
             %% make separate u and v from velocity vector u_crt
             %
             u_x = u_crt(:,1); u_y = u_crt(:,2);
             %
             %% write the flow-field data
             %
             filname = "2DFlowField";
             fnameip1 = sprintf(filname+".txt");
             fileID = fopen(fnameip1,'w');
             %
             AA = [u_x, u_y, p_crt];
             AA_new = AA';
             %
             fprintf(fileID,'%10.5e %10.5e %10.5e\n',AA_new);
             fclose(fileID);
             %
             %% write the time and cnt 
             %
             BB = [time ; cnt; p ; m];
             %
             filname = "timeANDcount";
             fnameip2 = sprintf(filname+".txt");
             fileID = fopen(fnameip2,'w');
             fprintf(fileID,'%10.5e %10.5d %5d %5d\n',BB);
             fclose(fileID);
             %
         end
         %
         disp('Hola')
        %% Local postprocess for superconvergence
        % (Warning: only straight-sided elements!)
        %{
        disp('Performing local postprocess...')
        p_star = referenceElement.degree + 1;
        nOfNodes_star = 0.5*(p_star+1)*(p_star+2);
        referenceElement_star = createReferenceElementTriSuperParametric(p_star,degreeGeo);
[u_star,shapeFunctionss] = HDGpostprocessStokes(X,T,u,q,referenceElement_star,referenceElement,Elements,muElem(1),example);
        disp('Done!')
        %}
       
        ErrorCalculations
        
        %{
        
        % Relative error for the postprocessed solution of u:
        error_postd1 = zeros(length(Elements.D1),1);
        coordRef_star = referenceElement_star.NodesCoord;
        
        for i = 1:length(Elements.D1)
            iElem=Elements.D1(i);
            Te_lin = T(iElem,:);
            Xold = X(Te_lin,:);
            Xe=shapeFunctions*Xold;
            ind = (iElem-1)*nOfNodes_star+1:iElem*nOfNodes_star;
            error_postd1(iElem) = computeL2Norm(referenceElement_star,Xe,(1:nOfNodes_star),u_star(ind),u0,time);
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
            error_postint(iElem) = computeL2Norm_cut(LSe_star,referenceElement_star,Xe,(1:nOfNodes_star),u_star(ind),u0,time);
        end
        ErrorPostInt = sqrt(sum(error_postint.^2));
        %disp(['Error HDG postprocessed Int= ', num2str(ErrorPostInt)]);

        
        
        
        Error_global_post=sqrt(ErrorPostD1^2+ErrorPostInt^2);
        disp(['Error HDG postprocessed u= ', num2str(Error_global_post)]);
        disp(' ')
        %}
        
        
        
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

       
        
        %Plot solution
        disp('Hola')
        figure(22),clf
        PlotDiscontSol(u_crt(:,1),X,T,LS,referenceElement,Elements)
        colorbar, title('X-HDG numerical solution of u1')
        %
        disp('Hola')
        %
        figure(33),clf
        PlotDiscontSol(u_crt(:,2),X,T,LS,referenceElement,Elements)
        colorbar, title('X-HDG numerical solution of u2')
        %
        disp('Hola')
        %        
        figure(44),clf
        PlotDiscontSol(p_crt,X,T,LS,referenceElement,Elements)
        colorbar, title('X-HDG numerical solution of pressure')
        %
        disp('Hola')
        %
        figure(55),clf
        plotContinuosSolution(X,T,analyticalpressure(X,time,Re),referenceElement)
        colorbar, title('HDG solution analytical of pressure')
        %caxis([-1 1]);
       
        
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
        errors = [errors, ErrorP_global];
        %errorsq = [errorsq, Errorq_global];
        %errorsPost = [errorsPost, Error_global_post]; 
        
    end
    %{
    %% Other plots
    figure(21), clf
    plot(log10(hs),log10(errors),'-o',log10(hs),log10(errorsPost),'o-');
    legend('u','u*')

    slopes = (log10(errors(2:end))-log10(errors(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
    %slopesq = (log10(errorsq(2:end))-log10(errorsq(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
    slopesPost = (log10(errorsPost(2:end))-log10(errorsPost(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
 

    errors1=[errors1 ; errors];    
    errorspost1=[errorspost1 ; errorsPost];
    lhs = log10(hs);
    %}
    
 end