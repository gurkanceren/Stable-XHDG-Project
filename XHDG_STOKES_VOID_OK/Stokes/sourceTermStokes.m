function s = sourceTermStokes(Xg,mu)

x=Xg(:,1); y=Xg(:,2);

%s=Xg*0+1; %linear solution p=x+y+C
%s=[0.*y , 0.*y]; %p=C
%s=[1-mu*6*x , 1+6*mu*y]; %cubic solution with p=x+y+C
%s=[1-mu.*100.*x.^3 , 1-300.*y.*x.^2.*-mu]; 
s=[1-mu.*20.*y.^3 , 1-mu.*20.*x.^3];
