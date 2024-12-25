function p0=DirichletConditionPressCorr(X,time,Re)
   p0 = analyticalpressure(X,time,Re);
   %x = X(:,1);
   %y = X(:,2);
   %u0 =(0.*x)+(0.*y);
   %if (x>=-0.6)&(x<=0.6)
   %    u0 = analiticalSolutionLaplace(X);
   %else
   %    u0 =(0.*x)+(0.*y);
   %end
   %disp('Hola');