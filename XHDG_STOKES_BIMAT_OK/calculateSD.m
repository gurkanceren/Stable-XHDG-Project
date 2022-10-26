function [sD1,sD2]=calculateSD(igauss)
global mu1 mu2

x=igauss(:,1);
y=igauss(:,2);

sD1=[cos(1-y^5)-y^5] ; 
sD2=[sin(1-x^5)-x^5] ;
