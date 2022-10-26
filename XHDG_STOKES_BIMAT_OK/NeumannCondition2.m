
function [ge]=NeumannCondition2(n_g,X)
global mu1 mu2 
%velocity component 2
x=X(:,1);
y=X(:,2);

gradux1=5.*x^4 ;  %1 and 2 at the end of the components here refer to domain numbers
gradux2=-5*x^4*cos(x^5-1);
graduy1=x.*0;
graduy2=x.*0;

ge=n_g*[-mu1*[gradux1;graduy1] + mu2*[gradux2;graduy2]]+n_g*[x.*0+0 ; x.*0+0];





