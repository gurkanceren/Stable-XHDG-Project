function LS = EvaluateLS(X, example)
% 
% LS = EvaluateLS(X, example)
% LS is the value of the level set function at points X


if example == 1
    P0 = [0,0];  
    r0 = 0.5; 
    LS = sqrt((X(:,1) - P0(1)).^2 + (X(:,2) - P0(2)).^2) - r0;
else
   
    
     %y=mx+b straight line 
    
%     m = -2.8;  
%     b = -1; 
%     LS = m*(X(:,1)) + b - X(:,2);

    LS=-(X(:,1)-0.2);
    
%     for i=1:length(X(:,1))
%         
%         if X(i,1)>=0.75 || X(i,1)<=-0.75
%             
%             LS(i)=abs(LS(i));
%         end
%     end
            
end


