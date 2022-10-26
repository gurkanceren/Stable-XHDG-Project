function [L2error,L2norm] = computeL2ErrorHDG(analyticalSolution,u,referenceElement,X,T)
%[L2error,L2norm] = computeL2ErrorHDG(@analyticalSolution,u,referenceElement,X,T)

[nOfElements,nOfElementNodes] = size(T);
%Information of the reference element
IPw = referenceElement.IPweights; 
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;
%Number of Gauss points
ngauss = length(IPw);

L2norm = 0;
L2error = 0;

%Loop in elements
for iElem = 1:nOfElements
    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    ue = u(ind,:);
    Xe = X(T(iElem,:),:); xe = Xe(:,1); ye=Xe(:,2);
    Xg = N*Xe;
    ua = feval(analyticalSolution,Xg);
    uh = N*ue;
    %Loop in Gauss points
    for g = 1:ngauss
        %Values at current integration point
        N_g = N(g,:); Nxi_g = Nxi(g,:); Neta_g = Neta(g,:);
        %Jacobian
        J = [Nxi_g*xe	  Nxi_g*ye
            Neta_g*xe  Neta_g*ye];
        %Integration weight
        dvolu=IPw(g)*det(J);
        %Contribution of the current integration point to the elemental L2 Norm
        L2error = L2error + sum((uh(g,:)-ua(g,:)).^2)*dvolu;
        L2norm = L2norm + sum((ua(g,:)).^2)*dvolu;
    end
end
L2error = sqrt(L2error);
L2norm = sqrt(L2norm);

