%% NumSimLab: 18-Feb-2024
%
%% Academic 2D HDG code for numerical simulation of 
%% the unsteady incompressible Navier-Stokes equations 
%% subject to Dirichlet boundary conditions.
%
%% *****************
%% This code uses 3rd order Backward Difference Scheme for time integration.
%% *****************
%
%% Modified from 2D HDG code for unsteady linear Conv-Diffn. eqn. with Dirchlet B.C.
%% Developed by Dr. Haroon Ahmad under the leadership of Dr. Ceren Gurkan. 
%% Dr. Ceren Gurkan (Asst. Prof.,Civil Engg. Deptt.,Kadir Has University,Istanbul,Turkey).
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
%% uses the HDG projection-Step method...
%% ... by Ueckermann & Lermusiaux (JCP,306,(2016),Pg:390-421).
%
%
clc, clear all, close all, %home
restoredefaultpath, setpath

% initialize time:
time = 0;
% time-integration count:
cnt = 0;
% time-step:
dt = 5.0e-5;
% Reynolds number:
Re = 1;
% parameter used in discretized temporal/unsteady term:
%a_parm = 1;
% Stabilization parameter for Predictor Step:
tau = 1/Re;
% final time:
T_final = 0.1; %pi/4;
%*****************************

errors1=[];   errorspost1=[];

for p=1
    %degree
    errors = []; errorsPost = []; errorsq= []; hs=[];
    for m=4
        %mesh_number
        
        filename = ['mesh' num2str(m) '_P' num2str(p) ];
        display('Solving...')
        fname1=[filename '.dcm'];
        display(fname1)
        
        hs=[hs,0.5^(m+1)];

        % Load data
        meshName = fname1;
        if all(meshName(end-2:end)=='dcm')
            GenerateMatFileFromEZ4U(['Meshes3/' meshName]);
        end
        load(['Meshes3/' meshName(1:end-3) 'mat']);

%X = (X + 1)./2; % modify the mesh to have it defined in [0,1] when using Meshes3
%X = (2*X-1);  % modify the mesh to have it defined in [-1,1]
%X = X+1; % modify the mesh to have it defined in [0,2] when using Meshes3
figure(1),clf
plotMesh(X,T)
figure(1), hold on; 
%AddFrontPlot (2,X);
%figure(1), hold on; 
%AddElementNumber(X,T) 
%Reference element
referenceElement = createReferenceElement(1,elemInfo.nOfNodes);

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
%p_wrc_crt = zeros(nOfElements*nOfElementNodes,1);

%allocate viscosity in the elements; 
%mu1=1;
%mu2=1;

%Viscosity 

%% Computation
time=0;
%% the time-loop starts here:
%while(cnt<100) 
while(time<T_final)

    time = time+dt; % update time
    cnt = cnt + 1;
    %
    % save previous time data
    %u_old = u_crt;
    p_old = p_crt;
    %
    %storage of u values of old times
    if(cnt==1)      %for n=1
        u_old = u_crt;  
        %here: u_old = u_(n-1)
        % multiplier for discretized unsteady velocity term at n time step
        a_parm = 1;     
    elseif(cnt==2)  %for n=2
        u_old2 = u_crt; 
        %here: u_old2 = u_(n-1) and u_old = u_(n-2)
        % multiplier for discretized unsteady velocity term at n time step
        a_parm = 2/3;       
    elseif(cnt==3)  %for n=3
        u_old3 = u_crt; 
        %here: u_old3 = u_(n-1) and u_old2 = u_(n-2) and u_old = u_(n-3)
        % multiplier for discretized unsteady velocity term at n time step
        a_parm = 6/11;        
    else            %for n>3
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
    %% Computation Step-1 (Predictor Step: gives Predictor velocity):
    % Loop in elementsnOfFaces
%    disp('Predictor Step: Loop in elements...')
    %
    [K f QQ UU Qf Uf uDirichlet] = hdg_matrix_predictorStep(X,T,F,referenceElement,infoFaces,tau,dt,cnt,u_old,u_old2,u_old3,p_old,time,Re,a_parm);
    %
    % Face solution
%    disp('Solving linear system for predictor velocity: u...')
    lambda = K\f;
    uhat = [lambda; uDirichlet];
    %
    % Elemental solution
%    disp('Calculating element by element solution for predictor velocity...')
    [u,q]=computeElementsSolutionPredictorStep(uhat,UU,QQ,Uf,Qf,T,F,referenceElement);
    %
    
    %% Computation Step-2 (Correction-Pressure Poisson eqn.: gives correction-pressure):
    % Loop in elementsnOfFaces
%    disp('Correction Pressure Step: Loop in elements...')
    %
    [Kp fp QQp UUp Qfp Ufp pDirichlet] = hdg_matrix_PressCorr(X,T,F,referenceElement,infoFaces,tauP,dt,Re,a_parm,u,uhat,time);
    %
    % Face solution
%    disp('Solving linear system for correction-pressure: pcorr...')
    lambda_pcorr = Kp\fp;
    phat = [lambda_pcorr; pDirichlet];
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
    [u_crt,p_crt,p_wrc_crt]=computeElementalCorrectionStep_Qcorr(u,uhat,pcorr,phat,qcorr,X,T,F,referenceElement,infoFaces,tauP,dt,Re,a_parm,time,cnt,p_old);
    
    %disp('count of time-integration')
    disp(cnt)
    %
end  
%% the time-loop ends here.
%
%disp('Hola')
Xelemental=zeros(nOfElementNodes*nOfElements,2);
for i=1:nOfElements
    %
    Te=T(i,:);
    Xe=X(Te,:);
    ind=(i-1)*(nOfElementNodes)+(1:nOfElementNodes);
    Xelemental(ind,:)=Xe;        
end
%
p_analy = analyticalpressure(Xelemental,time,Re);
delta_p = p_crt-p_analy;
%% plots HDG solution
figure(2),clf
%plotDiscontinuosSolutionSuperParametric(X,T,u(:,1),referenceElement,20)
plotDiscontinuosSolution(X,T,u_crt(:,1),referenceElement,20)
colorbar, title('HDG solution: u1 Navier-Stokes')

figure(3),clf
%plotDiscontinuosSolutionSuperParametric(X,T,u(:,2),referenceElement,20)
plotDiscontinuosSolution(X,T,u_crt(:,2),referenceElement,20)
colorbar, title('HDG solution: u2 Navier-Stokes')

figure(4),clf
%plotDiscontinuosSolutionSuperParametric(X,T,p,referenceElement,20)
plotDiscontinuosSolution(X,T,p_crt,referenceElement,20)
colorbar, title('HDG solution: p Navier-Stokes')

figure(5),clf
%plotDiscontinuosSolutionSuperParametric(X,T,p,referenceElement,20)
plotDiscontinuosSolution(X,T,pcorr,referenceElement,20)
colorbar, title('HDG solution: pCorr Navier-Stokes')

figure(6),clf
%plotDiscontinuosSolutionSuperParametric(X,T,p,referenceElement,20)
plotDiscontinuosSolution(X,T,p_analy,referenceElement,20)
colorbar, title('HDG solution: analytical p Navier-Stokes')

figure(7),clf
%plotDiscontinuosSolutionSuperParametric(X,T,p,referenceElement,20)
plotDiscontinuosSolution(X,T,delta_p,referenceElement,20)
colorbar, title('HDG solution: delta_P')


%
%
disp('Performing local postprocess on corrected velocity...')
p_star = referenceElement.degree + 1;
nOfNodes_star = 0.5*(p_star+1)*(p_star+2);
referenceElement_star = createReferenceElement(1,nOfNodes_star);
[u_star,shapeFunctions] = hdg_postprocess(X,T,u_crt,q,referenceElement_star,referenceElement,Re);
disp('Done!')
%

%% Plots postprocess solution
%figure(22),clf
%plotDiscontinuosSolution(X,T,u_star(:,1),referenceElement_star,20)
%colorbar, title('post-processed HDG solution: u1*')

%figure(32),clf
%plotDiscontinuosSolution(X,T,u_star(:,2),referenceElement_star,20)
%colorbar, title('post-processed HDG solution: u2*')

%% Error

% analytical solution
%% Relative error for corrected velocity:

%% FROM HERE NOW:******

u1_ex = @analyticalVelocityStokes_u1;
u2_ex = @analyticalVelocityStokes_u2;
p_ex = @analyticalpressure;

% Relative error: u
error = zeros(nOfElements,1);
for iElem = 1:nOfElements
    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    error(iElem) = computeL2Norm(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u_crt(ind,1),u1_ex,time,Re)^2 ...
    +computeL2Norm(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),u_crt(ind,2),u2_ex,time,Re)^2;
end
Error = sqrt(sum(error));
disp(['Error HDG = ', num2str(Error)]);

% Relative error for the postprocessed solution: u
error_post = zeros(nOfElements,1);
coordRef_star = referenceElement_star.NodesCoord;
for iElem = 1:nOfElements
    Te_lin = T(iElem,:);
    Xold = X(Te_lin,:);
    Xe=shapeFunctions*Xold;
    ind = (iElem-1)*nOfNodes_star+1:iElem*nOfNodes_star;
    error_post(iElem) = computeL2Norm(referenceElement_star,Xe,(1:nOfNodes_star),u_star(ind,1),u1_ex,time,Re).^2 ...
        +computeL2Norm(referenceElement_star,Xe,(1:nOfNodes_star),u_star(ind,2),u2_ex,time,Re).^2;
end
ErrorPost = sqrt(sum(error_post));
disp(['Error HDG postprocessed = ', num2str(ErrorPost)]);
disp(' ')

% Relative error: p
error_p = zeros(nOfElements,1);
for iElem = 1:nOfElements
    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    error_p(iElem) = computeL2Normscalar(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),p_crt(ind),p_ex,time,Re);
end
Error_p = sqrt(sum(error_p.^2));
disp(['Error HDG pressure = ', num2str(Error_p)]);

% Relative error: p w/o rotational correction
error_p = zeros(nOfElements,1);
for iElem = 1:nOfElements
    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    error_p(iElem) = computeL2Normscalar(referenceElement,X(T(iElem,:),:),(1:nOfElementNodes),p_wrc_crt(ind),p_ex,time,Re);
end
Error_p = sqrt(sum(error_p.^2));
disp(['Error HDG pressure w/o rot. corr. = ', num2str(Error_p)]);

%
%
 %% ConvergencePlots
 errors = [errors, Error];
 errorsPost = [errorsPost, ErrorPost]; 
        
    end

    %% Other plots
figure(20), clf
plot(log10(hs),log10(errors),'-o',log10(hs),log10(errorsPost),'o-');
legend('u','u*')

slopesu = (log10(errors(2:end))-log10(errors(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
slopesuPost = (log10(errorsPost(2:end))-log10(errorsPost(1:end-1)))./(log10(hs(2:end))-log10(hs(1:end-1)))
 
errors1=[errors1 ; errors];    
errorspost1=[errorspost1 ; errorsPost];
lhs = log10(hs);

end





