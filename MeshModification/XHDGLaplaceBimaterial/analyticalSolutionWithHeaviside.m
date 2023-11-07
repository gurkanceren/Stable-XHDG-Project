function sol=analyticalSolutionWithHeaviside(X)

x = X(:,1);
y = X(:,2);
r=sqrt(x.^2+y.^2);

sol = [3*r.^6+((2*0.4^6)/2);+1.0*r.^6-(((2*0.4^6))/2)];