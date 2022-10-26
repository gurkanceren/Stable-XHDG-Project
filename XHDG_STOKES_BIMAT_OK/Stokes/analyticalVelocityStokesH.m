function vel=analyticalVelocityStokesH(X)

x=X(:,1); y=X(:,2);
vel=zeros(x.*0,4);
global mu1 mu2;

for i=1:length(y)
    vel (i,:) = [  ((x(i).^2*y(i))*(mu1+mu2)+(mu1*mu2)) /(2*mu1*mu2) ; ...
        ((x(i).^2*y(i))*(mu2-mu1)-(mu1*mu2)) /(2*mu1*mu2);...
        ((-x(i)*y(i)^2)*(mu1+mu2)+(mu1*mu2))/(2*mu1*mu2) ;...
        ((-x(i)*y(i)^2)*(mu2-mu1)-(mu1*mu2))/(2*mu1*mu2)];
end


% for i=1:length(y)
%     vel (i,:) = [(x(i).^2*y(i))*(mu1+mu2)/(2*mu1*mu2) ; (x(i).^2*y(i))*(mu2-mu1)/(2*mu1*mu2);...
%         ((-x(i)*y(i)^2)*(mu1+mu2)+(2*mu1*mu2))/(2*mu1*mu2) ; (x(i)*y(i)^2)*(-mu2+mu1)/(2*mu1*mu2)];
% end


% for i=1:length(y)
%     
%     if y(i)<= 0
%         %vel = [(x.^2*y)/mu1 ; y.*0 ; -x.*(y.^2)/ ; y.*0]; %cubic solution
%         %vel (i) = [1-exp(Lmbd1)*sin(2*pi.*y(i))];  % ; y.*0; y.*0 ; y.*0];
%         vel(i,:) = [(x(i).^2*y(i))/mu1 ; y(i)*0 ; (-x(i)*y(i)^2)/mu1 ; y(i)*0]; %cubic solution
% 
%     else 
%         %vel (i) = [1-exp(Lmbd2)*sin(2*pi.*y(i))];  % ; y.*0; y.*0 ; y.*0];
%         %vel = [x.^3 ; y.*0 ; -3*(x.^2).*y ; y.*0]; %cubic solution
%         vel(i,:) = [(x(i).^2*y(i))/mu2 ; y(i)*0 ; (-x(i)*y(i)^2)/mu2 ; y(i)*0];
%     end
% end

%vel = [x.^3 ; y.*0 ; -3*(x.^2).*y ; y.*0]; %cubic solution

%vel = [x.^3 ; y.*0 ; -3*(x.^2).*y ; y.*0]; %cubic solution

%vel=[vel ; vel.*0 ; vel.*0 ; vel.*0];

%vel = [(x+y).*0+5 ;(x+y).*0; (x-y).*0+5; (x-y).*0] ;%(x+y)*0 ; (x-y)*0 ]; %linear solution
%vel = [((x+y).*0)+3 ;((x+y).*0)+1; ((x+y).*0)+5; ((x-y).*0)+1] ;%(x+y)*0 ; (x-y)*0 ]; %linear solution
%vel = [0.*y+5; 0.*x+5];
%vel = [x.^3 ; -3*(x.^2).*y]; %cubic solution
%vel = [y.^3 ; x.^3]; 
