

function [neux] = neumannu1 (nx,ny,I,mu) 


gradu=zeros (size(I(:,1),1),2);
gradu (:,1)=gradu (:,1)+0;
gradu (:,2)=gradu (:,2)+3.*I(:,2).^2;
p=[I(:,1)+I(:,2)]; %p=x+y+c
normal=[nx,ny];
for i=1:length(nx)   
   one(i,1)=-mu*gradu(i,:)*normal(i,:)';   
end
neux= one + p.*nx;