function  [u,q,u_analy]=computeElementsSolution(lambda,UU,QQ,Uf,Qf,T,F,X,time)
 
nOfElements = size(T,1);
nOfElementNodes = size(T,2);
nOfFaceNodes = size(UU{1},2)/3;

u = zeros(nOfElements*nOfElementNodes,1);
q = zeros(2*nOfElements*nOfElementNodes,1);
u_analy = zeros(nOfElements*nOfElementNodes,1);


%Loop in elements
for ielem= 1:nOfElements
    %
    Te = T(ielem,:);
    Xe = X(Te,:);
    %
    Fe = F(ielem,:);
    aux = (1:nOfFaceNodes);
    ind = [(Fe(1)-1)*nOfFaceNodes + aux,(Fe(2)-1)*nOfFaceNodes + aux,(Fe(3)-1)*nOfFaceNodes + aux];
    u((ielem-1)*nOfElementNodes+(1:nOfElementNodes)) = UU{ielem}*lambda(ind) + Uf{ielem};
    q((ielem-1)*2*nOfElementNodes+(1:2*nOfElementNodes)) = QQ{ielem}*lambda(ind) + Qf{ielem};
    %a = analiticalSolutionLaplace(Xe,time);
    %b = UU{ielem}*lambda(ind) + Uf{ielem};
    %c = (ielem-1)*nOfElementNodes+(1:nOfElementNodes);
    u_analy((ielem-1)*nOfElementNodes+(1:nOfElementNodes)) = analiticalSolutionLaplace(Xe,time);
    %disp('Hola')   
end


