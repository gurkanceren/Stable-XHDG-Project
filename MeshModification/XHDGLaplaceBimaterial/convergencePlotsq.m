


function convergencePlotsq

clc 

errors=[63.1576  30.2455  10.2992  3.0421
         23.4498  5.4821  0.93297  0.15287 
          5.6228  0.58498  0.048977  0.0040917 
             0.96567  0.03849  0.0015712   6.6461e-05];
    
    

  
hs = 0.5.^[1 2 3 4];
lerrors = log10(errors);
lhs = log10(hs);


figure
semilogy(lhs (1,:),errors(1,:),'-bo','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errors(2,:),'-rs','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errors(3,:),'-g<','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errors(4,:),'-md','LineWidth',3,'MarkerSize',12)
hold on

[h]=legend('P1','P2','P3','P4','Location','southeast');
%h=findobj(ax,'type','axis');
%set(h(1), 'visible', 'off')    % Hides the legend's axes (legend border and background)
set(h,'PlotBoxAspectRatioMode','manual');
set(h,'PlotBoxAspectRatio',[1 1.1 1]);
set(h,'FontSize',14)
%set(icons,'FontSize',30)
%set(str,'FontSize',30)
xlabel('log_{10}(Mesh Size)','FontSize',20,'fontname','times','fontweight','bold'), ylabel('L2 error Q','FontSize',20,'fontname','times','fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')
set(gca,'fontname','times')
%title('HDG','FontSize',15)
%legend boxoff 


p1slope1=(log10(errors(1,1))-log10(errors(1,2)))/(log10(0.5^1)-log10(0.5^2))
p1slope2=(log10(errors(1,2))-log10(errors(1,3)))/(log10(0.5^2)-log10(0.5^3))
p1slope3=(log10(errors(1,3))-log10(errors(1,4)))/(log10(0.5^3)-log10(0.5^4))

p2slope1=(log10(errors(2,1))-log10(errors(2,2)))/(log10(0.5^1)-log10(0.5^2))
p2slope2=(log10(errors(2,2))-log10(errors(2,3)))/(log10(0.5^2)-log10(0.5^3))
p2slope3=(log10(errors(2,3))-log10(errors(2,4)))/(log10(0.5^3)-log10(0.5^4))

p3slope1=(log10(errors(3,1))-log10(errors(3,2)))/(log10(0.5^1)-log10(0.5^2))
p3slope2=(log10(errors(3,2))-log10(errors(3,3)))/(log10(0.5^2)-log10(0.5^3))
p3slope3=(log10(errors(3,3))-log10(errors(3,4)))/(log10(0.5^3)-log10(0.5^4))

p4slope1=(log10(errors(4,1))-log10(errors(4,2)))/(log10(0.5^1)-log10(0.5^2))
p4slope2=(log10(errors(4,2))-log10(errors(4,3)))/(log10(0.5^2)-log10(0.5^3))
p4slope3=(log10(errors(4,3))-log10(errors(4,4)))/(log10(0.5^3)-log10(0.5^4))




