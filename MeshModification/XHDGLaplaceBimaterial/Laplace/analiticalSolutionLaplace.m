function u = analiticalSolutionLaplace(X)
global mu1 mu2 R
%%%%%% Peraire Example  

x = X(:,1);
y = X(:,2);

u = x*0;

for i=1:length(x)
    
    r(i)=sqrt(x(i)^2+y(i)^2);
    
    if r(i)>=R    %% Domain 1 
        
        u(i)=(1/mu1)*r(i)^5;
        
    elseif  r(i)<R
        
        u(i)=(1/mu2)*r(i)^5+((1/mu1)-(1/mu2))*R^5;
        
    end
    
end





