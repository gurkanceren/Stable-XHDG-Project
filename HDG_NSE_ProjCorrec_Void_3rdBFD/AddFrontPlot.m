
function AddFrontPlot (example,X)

Radius=0.5;

if example == 1
        P0 = [0,0];
        tt = linspace(0,pi*2);
        xx = P0(1) + Radius*cos(tt);
        yy = P0(2) + Radius*sin(tt);
        plot(xx,yy,'r--')

        
else
    
%     m = -2.8;  
%     b = -1;
%     
%     yy= m*(X(:,1)) + b;
%     xx= (yy-b)/m;
%     
%         plot(xx,yy,'r-')
 
yy=linspace(1,-1,100);
        
%yy=(line+0.2);

xx=ones(length(yy))*0.2;
plot(xx,yy,'r-','LineWidth',2)



end
