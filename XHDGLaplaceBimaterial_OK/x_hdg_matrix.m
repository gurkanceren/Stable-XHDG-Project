function [KK,f, QQ, UU, Qf, Uf] = x_hdg_matrix(X,T,F,referenceElement,infoFaces,tau,LS,Elements,mu1,mu2)
% This routine is valid only for triangles

nOfElements = size(T,1);
nOfElementNodes = size(T,2);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(infoFaces.intFaces,1)+size(infoFaces.extFaces,1);

nDOF = (nOfFaces)*2*nOfFaceNodes;  % as if all the faces were enriched !!!!
KK = spalloc(nDOF,nDOF,nOfFaceNodes*nDOF);
f = zeros(nDOF,1);

QQ = cell(nOfElements,1); UU = cell(nOfElements,1); Qf = cell(nOfElements,1); Uf = cell(nOfElements,1);

% %TEXT conservativity conditions (inicialization)
% residualFaces=zeros(2*nOfFaces*nOfFaceNodes,1); 
% allFaces = [infoFaces.intFaces(:,1:2);infoFaces.extFaces]; 
% gama = [];
% for iface=1:size(allFaces,1)
%     elem1 = allFaces(iface,1); faceElem1=allFaces(iface,2);
%     Tf = T(elem1,referenceElement.faceNodes(faceElem1,:)); Xf = X(Tf,:); LSf = LS(Tf);
%     if (all(LSf<0) | all(LSf>0)) %non-cut face
%         gamaF = [analiticalSolutionLaplace(Xf);zeros(nOfFaceNodes,1)];
%     else %cut face
%         gamaF = analyticalSolutionWithHeaviside(Xf);
%     end
%     gama = [gama;gamaF];
% end
%%end TEST inicialization

% loop in elements
for iElem = 1:nOfElements
    %disp(sprintf('\nELEMENT #%d: ',iElem))
    
    Te = T(iElem,:);
    Xe = X(Te,:);
    Fe = F(iElem,:);
    LSe = LS(Te);
    isFeInterior = (Fe <= nOfInteriorFaces); %Boolean (1=interior face)
    
    % elemental matrices
    d1=find(iElem==Elements.D1);    d2=find(iElem==Elements.D2);
    isCutElement = isempty([d1,d2]);
    if isCutElement %CUT element
        [Aqq,Auq,Aqu,Alu,Alq,fe,All,Auu]=xhdg_matrix_cut(referenceElement,...
            mu1,mu2,tau,Xe,LSe,iElem,F);
        nDOFface = 2*nOfFaceNodes; nDOFelement = 2*nOfElementNodes;
    else %standard element
        if isempty(d2), mu =mu1; else mu=mu2; end
        [Auq,Aqu,Alu,Alq,Auu,Aqq,All,fe]=xhdg_matrix(referenceElement...
            ,mu,tau,Xe,LSe,iElem,F);
        nDOFface = nOfFaceNodes; nDOFelement = nOfElementNodes;
    end
    % Local solver
    Aul = -Alu'; Aql = Alq';
    A = [Auu Auq; Aqu Aqq];
    UQ = -A\[Aul;Aql];
    fUQ= A\[fe;zeros(2*nDOFelement,1)];
    Ue = UQ(1:nDOFelement,:); Ufe=fUQ(1:nDOFelement); 
    Qe = UQ(nDOFelement+1:end,:); Qfe=fUQ(nDOFelement+1:end); 
    
    % Interior faces seen from the second element are flipped to have
    % proper orientation
    flipFace = zeros(1,3); %Boolean (1=face to be flipped)
    flipFace(isFeInterior)=any(infoFaces.intFaces(Fe(isFeInterior),[1 3])<iElem,2);
    aux=nOfFaceNodes:-1:1;
    if isCutElement
        aux = [aux,nOfFaceNodes+aux];
    end
    indL=1:3*nDOFface;
    indflip=[aux,nDOFface+aux,2*nDOFface+aux];
    aux=ones(1,nDOFface); aux = [flipFace(1)*aux, flipFace(2)*aux, flipFace(3)*aux];
    indL(aux==1)=indflip(aux==1); %permutation for local numbering
    
    Qe=Qe(:,indL);
    Ue=Ue(:,indL);
    Alq=Alq(indL,:);
    Alu=Alu(indL,:);
    All=All(indL,indL);
    
    %The local problem solver is stored for postprocess
    QQ{iElem} = sparse(Qe);  UU{iElem} = sparse(Ue);
    Qf{iElem} = sparse(Qfe); Uf{iElem} = sparse(Ufe);
    
    %Elemental matrices to be assembled
    KKe = Alq*Qe + Alu*Ue + 0.5*All;
    ffe = -(Alq*Qfe + Alu*Ufe);
    
    %All faces are assembled (Dirichlet included)
    aux = (1:nDOFface);
    ind = [(Fe(1)-1)*2*nOfFaceNodes + aux,(Fe(2)-1)*2*nOfFaceNodes + aux,(Fe(3)-1)*2*nOfFaceNodes + aux];       
    KK(ind,ind)=KK(ind,ind)+ KKe;
    f(ind) = f(ind) + ffe;
    
%     %TEST of conservativity condition at faces
%     if isCutElement
%         u_analy = analyticalSolutionWithHeaviside(Xe);
%     else
%         u_analy = analiticalSolutionLaplace(Xe);
%     end
%     q_analy = zeros(2*length(u_analy),1);
%     q_analy(1:2:2*nOfElementNodes)=-12*Xe(:,1).^2;
%     residualFaces(ind) = residualFaces(ind)+ Alu*u_analy + Alq*q_analy + 0.5*All*gama(ind);
%     %end TEST
%     
%     %TEST of local solver
%     ue = Ue*gama(ind)+Ufe; qe = Qe*gama(ind)+Qfe;
%     if max(abs([ue-u_analy;qe-q_analy]))>1.e-7
%         disp('Wrong local solver:')
%         u_and_uLS = [u_analy,ue], q_and_qLS = [q_analy,qe]
%     end
%     %end TEST

end
%figure(2), clf, spy(KK)

% %TEST output
% residualFaces = reshape(residualFaces,2*nOfFaceNodes,[]);
% maxAbsResidualInteriorFaces = max(abs(residualFaces(:,1:nOfInteriorFaces)))';
% rightFaces=find(maxAbsResidualInteriorFaces<1.e-10)
% facesWithWrongResiduals=find(maxAbsResidualInteriorFaces>1.e-10)
% elementsWrongFaces = infoFaces.intFaces(facesWithWrongResiduals,[1,3])
% %endTEST



        














