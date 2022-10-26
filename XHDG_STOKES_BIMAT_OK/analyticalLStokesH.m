function L_analy=analyticalLStokesH(X)

x=X(:,1); y=X(:,2);
L_analy=[];
global mu1 mu2;

for i=1:length(y)
    L_analy (i,:) = [(2*x(i)*y(i))*(mu1+mu2)/(2*mu1*mu2) ; (2*x(i)*y(i))*(mu2-mu1)/(2*mu1*mu2);...
        (x(i).^2)*(mu2+mu1)/(2*mu1*mu2)  ;  (x(i).^2)*(mu2-mu1)/(2*mu1*mu2);...
        (-y(i)^2)*(mu1+mu2)/(2*mu1*mu2) ; (y(i)^2)*(mu1-mu2)/(2*mu1*mu2);...
        (-2*x(i)*y(i))*(mu2+mu1)/(2*mu1*mu2)  ;  (2*x(i)*y(i))*(-mu2+mu1)/(2*mu1*mu2) ];
        
end
