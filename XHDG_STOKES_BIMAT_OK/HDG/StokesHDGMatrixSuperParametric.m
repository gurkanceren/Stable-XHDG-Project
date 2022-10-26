function [KK f localSolverMat localSolverVec] = StokesHDGMatrixSuperParametric(X,T,F,referenceElementSuperParametric,infoFaces,tau,Elements,LS)
% [KK f localSolverMat localSolverVec] = StokesHDGMatrixSuperParametric(muElem,X,T,F,referenceElementSuperParametric,infoFaces,tau)
% The geometry is approximated with the degree given by the mesh (X,T) and the
% reference element in referenceElementSuperParametric
% Nodal values are assumed to be ordered by components: u = [...u1... ...u2...]
global mu1 mu2;
nOfFaces = max(max(F));
nOfElements = size(T,1);
nOfElementNodes = size(T,2);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfFaceNodes = size(referenceElementSuperParametric.NodesCoord1d,1);
aux=nOfFaceNodes:-1:1; aux=[aux,aux+nOfFaceNodes]; 
indflip=[aux,2*nOfFaceNodes+aux,4*nOfFaceNodes+aux];
f = zeros(4*nOfFaces*nOfFaceNodes+nOfElements,1);
localSolverMat = cell(nOfElements,1);
localSolverVec = cell(nOfElements,1);
KK = sparse(nOfFaces*4*nOfFaceNodes+nOfElements,nOfFaces*4*nOfFaceNodes+nOfElements);
%KK=zeros(nOfFaces*4*nOfFaceNodes+nOfElements);

% loop in elements
nOfuhatNodes=4*nOfFaces*nOfFaceNodes;
%indK = 1; indf=1; n = nOfElements*(1+2*3*nOfFaceNodes)^2; ind_i  =
%zeros(1,n); ind_j  = zeros(1,n); coef_K = zeros(1,n);


% % %TEXT conservativity conditions (inicialization)
% % residualFaces=zeros(4*nOfFaces*nOfFaceNodes,1); 
% % allFaces = [infoFaces.intFaces(:,1:2);infoFaces.extFaces]; 
% % gama = [];
% % for iface=1:size(allFaces,1)
% %     elem1 = allFaces(iface,1); faceElem1=allFaces(iface,2);
% %     Tf = T(elem1,referenceElementSuperParametric.faceNodes(faceElem1,:)); Xf = X(Tf,:); LSf = LS(Tf);
% %     if (all(LSf<0) | all(LSf>0)) %non-cut face
% %         gamaF = [analyticalVelocityStokes(X)];
% %     else %cut face
% %         gamaF = analyticalVelocityStokesH(X);
% %     end
% %     gama = [gama;gamaF];
% % end
% % %end TEST inicialization
    
for iElem=1:nOfElements
    %display(iElem)
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    LSe= LS(Te);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)

    % elemental matrices
    d1=find(iElem==Elements.D1);    d2=find(iElem==Elements.D2);
    isCutElement = isempty([d1,d2]);
    
    if isCutElement %CUT element
       [localSolverMatElem,localSolverVecElem,Alu,AlL,Alp,All,Arl,AI] = KKeElementalMatricesSuperParametricCut...
        (Xe,referenceElementSuperParametric,tau(iElem,:),LSe);
         nDOFface = 4*nOfFaceNodes;
         aux=nOfFaceNodes:-1:1; aux=[aux,aux+nOfFaceNodes]; 
         indflip=[aux,2*nOfFaceNodes+aux,4*nOfFaceNodes+aux,6*nOfFaceNodes+aux,8*nOfFaceNodes+aux,10*nOfFaceNodes+aux];
         auxaux=[(1:nDOFface)];
         
    % Interior faces seen from the second element are flipped to have
    % proper orientation

    flipFace = zeros(1,3); %Boolean (1=face to be flipped)
    flipFace(isFeInterior)=any(infoFaces.intFaces(Fe(isFeInterior),[1 3])<iElem,2);
    indL=1:3*nDOFface; 
    aux=ones(1,nDOFface); aux = [flipFace(1)*aux, flipFace(2)*aux, flipFace(3)*aux];
    indL(aux==1)=indflip(aux==1); %permutation for local numbering
    indr = length(indL)+1;      
    localSolverMatElem=localSolverMatElem(:,[indL,indr]);
    AlL=AlL(indL,:);  Alu=Alu(indL,:); Alp=Alp(indL,:); All=All(indL,indL);
    Arl=Arl(:,indL);
        
    %The local problem solver is stored for postprocess
    localSolverMat{iElem} = sparse(localSolverMatElem);
    localSolverVec{iElem} = sparse(localSolverVecElem);
        
    %Elemental matrices to be assembled
    KKe = [Alu AlL Alp]*localSolverMatElem + [All,zeros(3*nDOFface,1)];
    KKe = [KKe ; [Arl, 0]];
    ffe = -[Alu AlL Alp]*localSolverVecElem;
    ffe = [ffe ; -AI];
         
         
    else  %standard element
        if isempty(d2), mu =mu1; else mu=mu2; end
    [localSolverMatElem,localSolverVecElem,Alu,AlL,Alp,All,Arl] = KKeElementalMatricesSuperParametric...
        (mu,Xe,referenceElementSuperParametric,tau(iElem,:)); 
        nDOFface = 2*nOfFaceNodes; 
         aux=nOfFaceNodes:-1:1; aux=[aux,aux+nOfFaceNodes]; 
         indflip=[aux,nDOFface+aux,2*nDOFface+aux];
         auxaux=[(1:nOfFaceNodes),(2*nOfFaceNodes+1:3*nOfFaceNodes)];
         
         
    % Interior faces seen from the second element are flipped to have
    % proper orientation

    flipFace = zeros(1,3); %Boolean (1=face to be flipped)
    flipFace(isFeInterior)=any(infoFaces.intFaces(Fe(isFeInterior),[1 3])<iElem,2);
    indL=1:3*nDOFface; 
    aux=ones(1,nDOFface); aux = [flipFace(1)*aux, flipFace(2)*aux, flipFace(3)*aux];
    indL(aux==1)=indflip(aux==1); %permutation for local numbering
    indr = length(indL)+1;      
    localSolverMatElem=localSolverMatElem(:,[indL,indr]);
    AlL=AlL(indL,:);  Alu=Alu(indL,:); Alp=Alp(indL,:); All=All(indL,indL);
    Arl=Arl(:,indL);
        
    %The local problem solver is stored for postprocess
    localSolverMat{iElem} = sparse(localSolverMatElem);
    localSolverVec{iElem} = sparse(localSolverVecElem);
        
    %Elemental matrices to be assembled
    KKe = [Alu AlL Alp]*localSolverMatElem + [All,zeros(3*nDOFface,1)];
    KKe = [KKe ; [Arl, 0]];
    ffe = -[Alu AlL Alp]*localSolverVecElem;
    ffe = [ffe ; 0];
         
         
    end    
    
    
    
    indrow = [(Fe(1)-1)*4*nOfFaceNodes + auxaux,(Fe(2)-1)*4*nOfFaceNodes + auxaux,(Fe(3)-1)*4*nOfFaceNodes + auxaux];
    indrow = [indrow, nOfuhatNodes+iElem]; %last component of Lambda_e=rho_e
    indcol = indrow;
    f(indrow) = f(indrow) + ffe;
    KK(indrow,indcol)=KK(indrow,indcol)+ KKe; 
    
        
%     %TEST of conservativity condition at faces
%     if isCutElement
%         u_analy = analyticalVelocityStokesH(Xe);
%     else
%         u_analy = analyticalVelocityStokes(Xe);
%     end
%     L_analy = zeros(2*length(u_analy),1);
%     p_analy = Xe(:,1)+ Xe(:,2);
%     residualFaces(ind) = residualFaces(ind)+ Alu*u_analy + AlL*L_analy + All*gama(ind) +Alp*p_analy;
    %end TEST
    
%     %TEST of local solver
%     display(iElem)
%     ue = U*gama(ind)+Uf; Le = Q*gama(ind)+Qf;
%     if max(abs([ue-u_analy;Le-L_analy]))>1.e-7
%         disp('Wrong local solver:')
%         u_and_uLS = [u_analy,ue], L_and_LLS = [L_analy,Le]
%     end
%     %end TEST
       
end

display('Hola')




