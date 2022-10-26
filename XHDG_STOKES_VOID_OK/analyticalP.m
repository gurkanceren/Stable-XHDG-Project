function vel=analyticalP(X,mu)

x=X(:,1); y=X(:,2);
%vel = [x+y ; x-y]; %linear solution
%vel = [0.*y+5; 0.*x+5];
%vel = [x.^3 ; -3*(x.^2).*y]; %cubic solution
vel = [x+y]; 
