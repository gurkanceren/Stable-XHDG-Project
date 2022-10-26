

function [u0] = AnalyticalSolution (X,T,Elements)

nOfElements=size(T,1);
nOfElementNodes=size(T,2);

u0 = zeros(nOfElements*2*nOfElementNodes,1);

for ielem=1:nOfElements
    
    d1=find(ielem==Elements.D1);    d2=find(ielem==Elements.D2);
    isCutElement = isempty([d1,d2]);
    Te=T(ielem,:);
    Xe=X(Te,:);
    Xf = [Xe;Xe];
    
    
    if ~isCutElement %standard element
        
        ind = (ielem-1)*(nOfElementNodes*2) + (1:2*nOfElementNodes);
        u0F = [analiticalSolutionLaplace(Xe); zeros(nOfElementNodes,1)];
        
    else  %cut element
        ind = (ielem-1)*(nOfElementNodes*2) + (1:2*nOfElementNodes);
        u0F = analyticalSolutionWithHeaviside(Xe);
        
    end
    u0(ind)=u0F;
end





