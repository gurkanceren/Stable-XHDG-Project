function u0=DirichletCondition(X)
%    u0 = analiticalSolutionLaplace(X);

u0=zeros(size(X,1),1);
for i=1:length(X)
    if X(i,1)<-0.2
     u0(i)= X(i,1)*2;
    else
     u0(i)= X(i,1)-0.2;
    end
    
end