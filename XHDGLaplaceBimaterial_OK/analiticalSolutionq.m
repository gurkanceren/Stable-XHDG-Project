function [qx, qy] = analiticalSolutionq(X)

% alpha = 0.1;
% beta = 0.3;
% a = 5.1;
% b = -6.2;
% c = 4.3;
% d = 3.4;

x = X(:,1);
y = X(:,2);

mu1=40;
mu2=1;

qx = x*0;
qy = y*0;

for i=1:length(x)
    
    if x(i)<0.1875
        
        qx(i)=6*(x(i)^5)*-mu1;
        qy(i)=0;
        
    else
        qx(i) = 240*(x(i)^5)*-mu2;
        qy(i)=0;
    end
    
end


