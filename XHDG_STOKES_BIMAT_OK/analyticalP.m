function pressure=analyticalP(X)

x=X(:,1); y=X(:,2);
pressure=x.*0;
% re1=1/mu1;
% re2=1/mu2;
% Lmbd1=(re1/2)-sqrt(((re1^2)/4)+4*pi^2);
% Lmbd2=(re2/2)-sqrt(((re2^2)/4)+4*pi^2);

pressure = [y+x ; x.*0]; 
 
% for i=1:length(y)    
%     if y<= 0.5
% 
%         pressure(i) = [0.5*exp(2*Lmbd1*x(i))]; 
%     else
% 
%         pressure(i) = [0.5*exp(2*Lmbd2*x(i))]; 
%     end
% end
% 
% 
%  pressure=[pressure ; pressure.*0];
