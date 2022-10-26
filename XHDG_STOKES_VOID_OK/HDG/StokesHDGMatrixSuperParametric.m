function [KK f localSolverMat localSolverVec] = StokesHDGMatrixSuperParametric(muElem,X,T,F,referenceElementSuperParametric,infoFaces,tau,Elements,LS)
% [KK f localSolverMat localSolverVec] = StokesHDGMatrixSuperParametric(muElem,X,T,F,referenceElementSuperParametric,infoFaces,tau)
% The geometry is approximated with the degree given by the mesh (X,T) and the
% reference element in referenceElementSuperParametric
% Nodal values are assumed to be ordered by components: u = [...u1... ...u2...]

nOfFaces = max(max(F));
nOfElements = size(T,1);
nOfElementNodes = size(T,2);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfFaceNodes = size(referenceElementSuperParametric.NodesCoord1d,1);
aux=nOfFaceNodes:-1:1; aux=[aux,aux+nOfFaceNodes]; 
indflip=[aux,2*nOfFaceNodes+aux,4*nOfFaceNodes+aux];
f = zeros(2*nOfFaces*nOfFaceNodes+nOfElements,1);
localSolverMat = cell(nOfElements,1);
localSolverVec = cell(nOfElements,1);
KK=sparse(nOfFaces*2*nOfFaceNodes+nOfElements,nOfFaces*2*nOfFaceNodes+nOfElements);

% loop in elements
nOfuhatNodes=2*nOfFaces*nOfFaceNodes;
%indK = 1; indf=1; n = nOfElements*(1+2*3*nOfFaceNodes)^2; ind_i  = zeros(1,n); ind_j  = zeros(1,n); coef_K = zeros(1,n);

for iElem=1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    LSe= LS(Te);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)

    % elemental matrices
    d1=find(iElem==Elements.D1);    d2=find(iElem==Elements.D2);
    isCutElement = isempty([d1]);
    if ~isempty(d2); continue; end;
    
    if isCutElement %CUT element
       [localSolverMatElem,localSolverVecElem,Alu,AlL,Alp,All,Arl] = KKeElementalMatricesSuperParametricCut...
        (muElem(1),Xe,referenceElementSuperParametric,tau(iElem,:),LSe);
        localSolverMatElem = [localSolverMatElem,[zeros(size(localSolverMatElem,1),1)]] ;

   elseif ~isempty(d1) %standard element
        
    [localSolverMatElem,localSolverVecElem,Alu,AlL,Alp,All,Arl] = KKeElementalMatricesSuperParametric...
        (muElem(iElem),Xe,referenceElementSuperParametric,tau(iElem,:));  
       
    end    
    
    % Interior faces seen from the second element are flipped to have
    % proper orientation
    flipFace = zeros(1,3); %Boolean (1=face to be flipped)
    flipFace(isFeInterior)=any(infoFaces.intFaces(Fe(isFeInterior),[1 3])<iElem,2);
    indL=1:6*nOfFaceNodes; 
    aux=ones(1,2*nOfFaceNodes); aux = [flipFace(1)*aux, flipFace(2)*aux, flipFace(3)*aux];
    indL(aux==1)=indflip(aux==1); %permutation for local numbering
    indr = length(indL)+1;      
    localSolverMatElem=localSolverMatElem(:,[indL,indr]);
    AlL=AlL(indL,:);  Alu=Alu(indL,:); Alp=Alp(indL,:); All=All(indL,indL);
    Arl=Arl(:,indL);
        
    %The local problem solver is stored for postprocess
    localSolverMat{iElem} = sparse(localSolverMatElem);
    localSolverVec{iElem} = sparse(localSolverVecElem);
        
    %Elemental matrices to be assembled
    KKe = [Alu AlL Alp]*localSolverMatElem + [All,zeros(2*3*nOfFaceNodes,1)];
    KKe = [KKe ; [Arl, 0]];
    ffe = -[Alu AlL Alp]*localSolverVecElem;
    ffe = [ffe ; 0];
    
    aux = (1:2*nOfFaceNodes);
    indrow = [(Fe(1)-1)*2*nOfFaceNodes + aux,(Fe(2)-1)*2*nOfFaceNodes + aux,(Fe(3)-1)*2*nOfFaceNodes + aux];
    indrow = [indrow, nOfuhatNodes+iElem]; %last component of Lambda_e=rho_e
    indcol = indrow;
    f(indrow) = f(indrow) + ffe;
    KK(indrow,indcol)=KK(indrow,indcol)+ KKe; 

       
end

for i=1:length(Elements.D2)
    ielem = Elements.D2(i);
localSolverMatD2= zeros(7*nOfElementNodes,3*2*nOfFaceNodes+1) ;
localSolverVecD2= zeros(7*nOfElementNodes,1) ;

localSolverMat{ielem} = sparse(localSolverMatD2);
localSolverVec{ielem} = sparse(localSolverVecD2);
end

