function [u,q]=computeElementsSolutionPredictorStep(lambda,UU,QQ,Uf,Qf,T,F,referenceElement)
%(uhat,rho,localSolverMat,localSolverVec,F)
%nOfElements = size(localSolverVec,1);
%[n,m]=size(localSolverMat{1});
%nOfElementNodes=n/7; nOfFaceNodes=(m-1)/6;
%
nOfElements = size(T,1);
nOfElementNodes = size(T,2);
%nOfFaceNodes = size(UU{1},2)/3;
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfFaces = size(referenceElement.faceNodes,1);

u = zeros(nOfElements*nOfElementNodes,2);
q = zeros(nOfElements*nOfElementNodes,4);
%u_Hat = zeros(nOfElements*nOfFaces*nOfFaceNodes,2);
%p = zeros(nOfElements*nOfElementNodes,1);

%Loop in elements
for ielem = 1:nOfElements
    Fe = F(ielem,:);
    aux = (1:2*nOfFaceNodes);
    ind = [(Fe(1)-1)*2*nOfFaceNodes + aux,(Fe(2)-1)*2*nOfFaceNodes + aux,(Fe(3)-1)*2*nOfFaceNodes + aux];
    %lambda = [uhat(ind);rho(ielem)];
    %sol = localSolverMat{ielem}*lambda + localSolverVec{ielem};    
    %sol = reshape(sol,nOfElementNodes,7);
    u_hat = lambda(ind);
    %
    sol_u = UU{ielem}*u_hat + Uf{ielem};
    %
    sol_u = reshape(sol_u,nOfElementNodes,2);
    %
    sol_q = QQ{ielem}*u_hat + Qf{ielem};
    %
    sol_q = reshape(sol_q,nOfElementNodes,4);
    %
    inf = (ielem-1)*nOfElementNodes+(1:nOfElementNodes);
    %
    u(inf,1) = sol_u(:,1); 
    u(inf,2) = sol_u(:,2);
    %
    q(inf,1) = sol_q(:,1); 
    q(inf,2) = sol_q(:,2); 
    q(inf,3) = sol_q(:,3); 
    q(inf,4) = sol_q(:,4);
    %
    %p(ind)=sol(:,7);
end
    %
    %disp('Hola')
    %

