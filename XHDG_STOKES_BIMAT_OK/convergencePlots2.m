
function convergencePlots2

%% Stokes convergence discnt  

errors=[ 0.9040    0.3892    0.1975    0.0536
        0.2811  0.1398   0.0051  6.4508e-04
         0.0217   0.0024  1.5592e-04   9.9888e-06 ] ;



errorspost=[ 0.2813    0.1392    0.0058    0.0008
           0.0341   0.0035   2.3644e-04   1.5557e-05
             0.0038  1.5763e-04  5.1704e-06   2.0426e-07 ];

         
         
%hs = 2./[2 4 8 16];
hs = 0.5.^[0 1 2 3];
lerrors = log10(errors);
lerrorspost = log10(errorspost);
lhs = log10(hs);

figure
x0=10;
y0=10;
width=500;
height=400;
set(gcf,'units','points','position',[x0,y0,width,height])
ylim([1e-8 1])


semilogy(lhs (1,:),errors(1,:),'-ko','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost(1,:),'--ko','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12)
hold on


semilogy(lhs (1,:),errors(2,:),'-ks','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost(2,:),'--ks','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12)
hold on


semilogy(lhs (1,:),errors(3,:),'-k<','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost(3,:),'--k<','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12)
hold on


% semilogy(lhs (1,:),errors(4,:),'-md','LineWidth',3,'MarkerSize',12)
% hold on
% semilogy(lhs (1,:),errorspost(4,:),'--md','LineWidth',3,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',12)
% hold on

[h]=legend('P1','P1 postp.','P2','P2 postp.','P3','P3 postp.','Location','southeast');
%[h]=legend('P2','P2 postp.','P3','P3 postp.','P4','P4 postp.','Location','southeast');
%h=findobj(ax,'type','axis');
%set(h(1), 'visible', 'off')    % Hides the legend's axes (legend border and background)
%set(h,'PlotBoxAspectRatioMode','manual');
%set(h,'PlotBoxAspectRatio',[1 0.8 1]);
set(h,'FontSize',20)
%set(icons,'FontSize',30)
%set(str,'FontSize',30)
xlabel('log_{10}(Mesh Size)','FontSize',20,'fontname','times','fontweight','bold'), ylabel('L2 error','FontSize',20,'fontname','times','fontweight','bold')
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

% p4slope1=(log10(errors(4,1))-log10(errors(4,2)))/(log10(0.5^1)-log10(0.5^2))
% p4slope2=(log10(errors(4,2))-log10(errors(4,3)))/(log10(0.5^2)-log10(0.5^3))
% p4slope3=(log10(errors(4,3))-log10(errors(4,4)))/(log10(0.5^3)-log10(0.5^4))


p1postslope1=(log10(errorspost(1,1))-log10(errorspost(1,2)))/(log10(0.5^1)-log10(0.5^2))
p1postslope2=(log10(errorspost(1,2))-log10(errorspost(1,3)))/(log10(0.5^2)-log10(0.5^3))
p1postslope3=(log10(errorspost(1,3))-log10(errorspost(1,4)))/(log10(0.5^3)-log10(0.5^4))

p2postslope1=(log10(errorspost(2,1))-log10(errorspost(2,2)))/(log10(0.5^1)-log10(0.5^2))
p2postslope2=(log10(errorspost(2,2))-log10(errorspost(2,3)))/(log10(0.5^2)-log10(0.5^3))
p2postslope3=(log10(errorspost(2,3))-log10(errorspost(2,4)))/(log10(0.5^3)-log10(0.5^4))

p3postslope1=(log10(errorspost(3,1))-log10(errorspost(3,2)))/(log10(0.5^1)-log10(0.5^2))
p3postslope2=(log10(errorspost(3,2))-log10(errorspost(3,3)))/(log10(0.5^2)-log10(0.5^3))
p3postslope3=(log10(errorspost(3,3))-log10(errorspost(3,4)))/(log10(0.5^3)-log10(0.5^4))

% p4postslope1=(log10(errorspost(4,1))-log10(errorspost(4,2)))/(log10(0.5^1)-log10(0.5^2))
% p4postslope2=(log10(errorspost(4,2))-log10(errorspost(4,3)))/(log10(0.5^2)-log10(0.5^3))
% p4postslope3=(log10(errorspost(4,3))-log10(errorspost(4,4)))/(log10(0.5^3)-log10(0.5^4))


% % hs = 0.5.^[3 4 5];
% % lhs = log10(hs);
% % figure (20)
% % 
% % semilogy(lhs (1,:),errors(1,:),'-ko','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errors(2,:),'-k>','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errors(3,:),'-ks','LineWidth',1)
% % 
% % semilogy(lhs (1,:),errorspost(1,:),'--ko','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errorspost(2,:),'--k>','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errorspost(3,:),'--ks','LineWidth',1)
% % 
% % legend('P2','P3','P4','P2 postprocesed','P3 postprocesed','P4 postprocesed')
% % xlabel('log_{10}(h)(Mesh Size)'), ylabel('L2 error')
% % title('Convergence Adapted HDG')




% % slopesFirstSegment = (lerrorspost(end,:)-lerrorspost(end-1,:))/(lhs(end)-lhs(end-1));
% % slopesSecondSegment = (lerrorspost(end-1,:)-lerrorspost(end-2,:))/(lhs(end-1)-lhs(end-2));
% %  SlopeTriangles (lhs,lerrorspost,slopesFirstSegment);
% %  
 