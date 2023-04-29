
function AddFrontPlot (example)

%Radius=0.21
Radius=0.41;

if example == 1
        %P0 = [0.5,0.5];
        P0 = [1,1];
        tt = linspace(0,pi*2);
        xx = P0(1) + Radius*cos(tt);
        yy = P0(2) + Radius*sin(tt);
        plot(xx,yy,'-','Color',[0 0 0],'LineWidth',3)
end

end
