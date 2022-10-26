function vel=analyticalVelocityStokes(X)
global mu1 mu2 R;
x=X(:,1); y=X(:,2);

for i=1:length(y)
    r(i)=sqrt(x(i)^2+y(i)^2);
    if r(i) >= R  %D1
        vel(i,:) = [y(i)^5 ; x(i)^5]; %cubic solution
    else 
        vel(i,:) = [cos(1-y(i)^5) ; sin(1-x(i)^5)];
    end
end


%vel = [(x+y).*0+5 ;(x+y).*0; (x-y).*0+5; (x-y).*0] ;%(x+y)*0 ; (x-y)*0 ]; %linear solution
%vel = [((x+y).*0)+3 ;((x+y).*0)+1; ((x+y).*0)+5; ((x-y).*0)+1] ;%(x+y)*0 ; (x-y)*0 ]; %linear solution
%vel = [0.*y+5; 0.*x+5];
%vel = [x.^3 ; -3*(x.^2).*y]; %cubic solution
%vel = [y.^3 ; x.^3]; 
