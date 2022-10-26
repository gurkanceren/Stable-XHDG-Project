function s = sourceTermStokes(Xg)
global mu1 mu2 R;
s=Xg.*0;
x=Xg(:,1); y=Xg(:,2);

for i=1:length(y)
   r(i)=sqrt(x(i)^2+y(i)^2); 
    if r(i)>= R  %D1

        s(i,:)=[1-mu1*20*y(i)^3 , 1-mu1*20*x(i)^3]; % with p=x+y+C with p=x^2+y^2+C
        
    else   

        %s(i,:)=[2*x(i)-mu1*20*y(i)^3 , 2*y(i)-mu1*20*x(i)^3];
        s(i,:)=[1-mu2*((-20*y(i)^3*sin(y(i)^5-1))-25*y(i)^8*cos(y(i)^5-1))  , 1-mu2*(25*x(i)^8*sin(x(i)^5-1)-20*x(i)^3*cos(x(i)^5-1))]; % with p=x+y+C
        
    end
end

