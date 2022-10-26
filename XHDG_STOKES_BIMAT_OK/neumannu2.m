

function [neuy] = neumannu2 (nx,ny,I,mu)

gradu=zeros (size (I(:,1),1),2);
gradu (:,1)=gradu (:,1)+3.*I(:,1).^2;
gradu (:,2)=gradu (:,2)+0;
p= [I(:,1)+I(:,2)]; %p=x+y+c
normal=[nx,ny];
for i=1:length(nx)    
   one(i,1)=-mu*gradu(i,:)*normal(i,:)';  
end
neuy= one + p.*ny;