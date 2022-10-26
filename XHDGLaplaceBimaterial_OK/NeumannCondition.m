
function [ge]=NeumannCondition(n_g,mu)


% delu/deln=gradu*n

gradu=[1;0];

ge=n_g*gradu*mu;





