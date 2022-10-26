function [ge]=NeumannCondition1(x,y)

%circle eqn--->f=sqrt(x^2+y^2-0.5)
%inner normal to integration points on I

if x==0 && y==0;
    
    ge=0;
    
else

ge=-y/sqrt(x.^2+y.^2);

end
