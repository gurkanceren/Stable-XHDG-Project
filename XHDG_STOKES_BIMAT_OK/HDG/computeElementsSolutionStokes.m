function [u,L,p]=computeElementsSolutionStokes(uhat,rho,localSolverMat,localSolverVec,F,Elements)
nOfElements = size(localSolverVec,1);
[n,m]=size(localSolverMat{Elements.Int(1)});
nOfElementNodes=n/7; nOfFaceNodes=(m-1)/6; %all enriched assumption

u = zeros(nOfElements*nOfElementNodes,2);
L = zeros(nOfElements*nOfElementNodes,4);
p = zeros(nOfElements*nOfElementNodes,1);

%Loop in elements
% for ielem = 1:nOfElements
%     Fe = F(ielem,:);
%     aux = (1:4*nOfFaceNodes);
%     ind = [(Fe(1)-1)*2*nOfFaceNodes + aux,(Fe(2)-1)*2*nOfFaceNodes + aux,(Fe(3)-1)*2*nOfFaceNodes + aux];
%     Lambda = [uhat(ind);rho(ielem)];
%     sol = localSolverMat{ielem}*Lambda + localSolverVec{ielem};
%     sol = reshape(sol,nOfElementNodes,7);
%     ind = (ielem-1)*nOfElementNodes+(1:nOfElementNodes);
%     u(ind,1) = sol(:,1); u(ind,2) = sol(:,2);
%     L(ind,1) = sol(:,3); L(ind,2) = sol(:,4); L(ind,3) = sol(:,5); L(ind,4) = sol(:,6);
%     p(ind)=sol(:,7);
% end

d1=[];
d2=[];

%Loop in elements
for ielem=1:nOfElements
    d1=find(ielem==Elements.D1);
    d2=find(ielem==Elements.D2);
    Fe = F(ielem,:);
    aux = (1:2*nOfFaceNodes);
    ind = [(Fe(1)-1)*2*nOfFaceNodes + aux,(Fe(2)-1)*2*nOfFaceNodes + aux,(Fe(3)-1)*2*nOfFaceNodes + aux]; 
    auxaux=[(1:(nOfFaceNodes/2)),(2*(nOfFaceNodes/2)+1:3*(nOfFaceNodes/2)), (4*(nOfFaceNodes/2)+1:5*(nOfFaceNodes/2)),...
        (6*(nOfFaceNodes/2)+1:7*(nOfFaceNodes/2)),(8*(nOfFaceNodes/2)+1:9*(nOfFaceNodes/2)), (10*(nOfFaceNodes/2)+1:11*(nOfFaceNodes/2))];
    indstd=ind(auxaux);
    if isempty(d1) && isempty(d2)  % cut element
    Lambda = [uhat(ind);rho(ielem)];
    sol = localSolverMat{ielem}*Lambda + localSolverVec{ielem};
    sol = reshape(sol,nOfElementNodes,7);
    ind = (ielem-1)*nOfElementNodes+(1:nOfElementNodes);
    u(ind,1) = sol(:,1); u(ind,2) = sol(:,2);
    L(ind,1) = sol(:,3); L(ind,2) = sol(:,4); L(ind,3) = sol(:,5); L(ind,4) = sol(:,6);
    p(ind)=sol(:,7);
    
    elseif ~isempty(d1) || ~isempty(d2)  %standard elements
    Lambda = [uhat(indstd);rho(ielem)]; 
    sol = localSolverMat{ielem}*Lambda + localSolverVec{ielem};
    sol = reshape(sol,(nOfElementNodes/2),7);
    ind = (ielem-1)*(nOfElementNodes)+(1:(nOfElementNodes/2));
    u(ind,1) = sol(:,1); u(ind,2) = sol(:,2);
    L(ind,1) = sol(:,3); L(ind,2) = sol(:,4); L(ind,3) = sol(:,5); L(ind,4) = sol(:,6);
    p(ind)=sol(:,7); 
    end
    
    d1=[];
    d2=[];
       
end