
function [ge]=NeumannCondition1(n_g,X)
global mu1 mu2 
%velocity component 1 
x=X(:,1);
y=X(:,2);


gradux1=x.*0;
gradux2=x.*0;
graduy1=5.*y.^4;
graduy2=-5*y^4*sin(y^5-1);

ge=n_g*[-mu1*[gradux1;graduy1] + mu2*[gradux2;graduy2]]+n_g*[x.*0+0;x.*0+0];





