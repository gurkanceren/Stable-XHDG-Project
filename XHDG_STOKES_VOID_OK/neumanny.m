

function [neuy] = neumannu2 (ny,I)

gradu=[0];   %change with different u 
p=I(:,1)+I(:,2);  % p=x+y+c
neuy= gradu.*ny + p.*ny;