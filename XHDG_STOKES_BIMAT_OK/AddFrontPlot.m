
function AddFrontPlot (example,X)

Radius=0.41;

if example == 1
        P0 = [0,0];
        tt = linspace(0,pi*2);
        xx = P0(1) + Radius*cos(tt);
        yy = P0(2) + Radius*sin(tt);
        plot(xx,yy,'k-','LineWidth',2)

        
else
    
%     m = -2.8;  
%     b = -1;
%     
%     yy= m*(X(:,1)) + b;
%     xx= (yy-b)/m;
%     
%         plot(xx,yy,'r-')
 
xx=linspace(-0.4,1.6,100);
        
%yy=(line+0.2);

yy=ones(length(xx))*0;
plot(xx,yy,'r-','LineWidth',2)



end
