
function AddFrontPlot (example,X)

Radius=0.5;

if example == 1
        P0 = [0,0];
        tt = linspace(0,pi*2);
        xx = P0(1) + Radius*cos(tt);
        yy = P0(2) + Radius*sin(tt);
        plot(xx,yy,'k-','LineWidth',2)

        
elseif example == 2
    
%     m = -2.8;  
%     b = -1;
%     
%     yy= m*(X(:,1)) + b;
%     xx= (yy-b)/m;
%     
%         plot(xx,yy,'r-')
 
yy=linspace(0,1,100);
        
%yy=(line+0.2);

xx=ones(length(yy))*0.85;
plot(xx,yy,'r-','LineWidth',2)



else
    
        [x,y] = meshgrid(-1:.02:1);
        R=sqrt(x.^2+y.^2);
        z=R-0.5-((y.^5 + 5*x.^4*y - 10*x.^2*y.^3)./(3*R.^4));
        surf(x,y,z)

end
