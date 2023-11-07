

function [ge]=NeumannConditiontest1(Pphy)

%line eqn--->
%inner normal to integration points on I

dy=Pphy(2,2)-Pphy(1,2);
dx=Pphy(2,1)-Pphy(1,1);
deltaux=-dy;
deltauy=dx;

ge=[deltauy];

end