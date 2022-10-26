function lap = analiticalLaplacianLaplace(X)
global R mu1 mu2 
%%%%% peraire Example

x = X(:,1);
y = X(:,2);
lap=0*x;

for i=1:length(lap)
    
r=sqrt(x(i)^2+y(i)^2);

if r>R % Domain 1
 
    lap(i)=10*r^3+15*x(i)^2*r+15*y(i)^2*r;
    %lap(i)=2;
else
    
   lap(i)=(1/mu2)*10*r^3+(1/mu2)*15*x(i)^2*r+(1/mu2)*15*y(i)^2*r;
   %lap(i)=2; 
    
end

end


    
    
    
    
    
    