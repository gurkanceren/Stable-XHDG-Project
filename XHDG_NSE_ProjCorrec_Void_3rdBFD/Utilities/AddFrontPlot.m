
function AddFrontPlot (example)

%Radius=0.25;
Radius=0.5;

if example == 1
        P0 = [1,0.5];
        %P0 = [0,0];
        tt = linspace(0,pi*2);
        xx = P0(1) + Radius*cos(tt);
        yy = P0(2) + Radius*sin(tt);
        plot(xx,yy,'-','Color',[0 0 0],'LineWidth',3)
        %
else
        %
        r0=0.37;
        r1=0.17;   
        [x,y]=meshgrid(-1:.01:1);
        theta = 2*atan2(y,x);
        A=(x).^2 + (y).^2 ;
        z=(sqrt(A) - r0 -r1.*cos(theta));
        contour(x,y,z,[0,-2],'-','Color',[0 0 0],'LineWidth',3)
        %
end

end
