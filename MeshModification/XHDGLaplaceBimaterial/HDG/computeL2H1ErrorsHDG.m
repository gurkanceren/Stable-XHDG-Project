function [L2error,H1error,L2norm,H1norm] = computeL2H1ErrorsHDG(analyticalSolutionAndDerivative,u,referenceElement,X,T)
%[L2error,H1error,L2norm,H1norm] = computeL2H1ErrorsHDG(@analyticalSolutionAndDerivative,u,referenceElement,X,T)
% u is assumed to be scalar

[nOfElements,nOfElementNodes] = size(T);
%Information of the reference element
IPw = referenceElement.IPweights; 
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;
%Number of Gauss points
ngauss = length(IPw);

L2norm = 0;
H1norm = 0;
L2error = 0;
H1error = 0;

%Loop in elements
for iElem = 1:nOfElements
    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    ue = u(ind);
    Xe = X(T(iElem,:),:); xe = Xe(:,1); ye=Xe(:,2);
    Xg = N*Xe;
    [ua,gradua] = feval(analyticalSolutionAndDerivative,Xg);
    uh = N*ue;
    duh_dxieta = [Nxi*ue Neta*ue];
    %Loop in Gauss points
    for g = 1:ngauss
        %Values at current integration point
        N_g = N(g,:); Nxi_g = Nxi(g,:); Neta_g = Neta(g,:);
        %Jacobian
        J = [Nxi_g*xe	  Nxi_g*ye
            Neta_g*xe  Neta_g*ye];
        %(x,y)-derivatives
        graduh_g = J\(duh_dxieta(g,:)'); 
        %Integration weight
        dvolu=IPw(g)*det(J);
        %Contribution of the current integration point to the elemental L2 Norm
        L2error = L2error + (uh(g)-ua(g))^2*dvolu;
        L2norm = L2norm + (ua(g))^2*dvolu;
        H1error = H1error + sum((graduh_g-gradua(g,:)').^2)*dvolu;
        H1norm = H1norm + sum(gradua(g,:).^2)*dvolu;
    end
end
H1error = sqrt(H1error+L2error);
L2error = sqrt(L2error);
H1norm = sqrt(L2norm+H1norm);
L2norm = sqrt(L2norm);


