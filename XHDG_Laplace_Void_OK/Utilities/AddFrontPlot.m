
function AddFrontPlot (example)

Radius=0.41;

if example == 1
        P0 = [0,0];
        tt = linspace(0,pi*2,50);
        xx = P0(1) + Radius*cos(tt);
        yy = P0(2) + Radius*sin(tt);
        plot(xx,yy,'-r','LineWidth',1,'Markersize',4)
end

end
