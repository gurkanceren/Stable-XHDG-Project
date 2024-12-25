function p0=DirichletConditionPressCorr(X)
   p0 = analiticalSolutionLaplace(X);
   %x = X(:,1);
   %y = X(:,2);
   %u0 =(0.*x)+(0.*y);
   %if (x>=-0.6)&(x<=0.6)
   %    u0 = analiticalSolutionLaplace(X);
   %else
   %    u0 =(0.*x)+(0.*y);
   %end
   %disp('Hola');