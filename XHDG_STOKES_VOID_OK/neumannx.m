

function [neux] = neumannu1 (nx,I)

gradu=[2];   %change with different u 
p=[I(:,1)+I(:,2)]; %p=x+y+c
neux= gradu.*nx + p.*nx;