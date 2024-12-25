
 mu_vector=zeros(nOfElements,1);

for iElem = 1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    d1=find(any(Xe(:,1)<0.2)==1);    
    d2=find(any(Xe(:,1)>0.2)==1);
    if isempty(d2)
    mu_vector(iElem)=mu1;
    else
    mu_vector(iElem)=mu2;    
    end
    
end