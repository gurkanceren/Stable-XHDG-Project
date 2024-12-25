

function convergencePlots

% % errors = [0.017186	0.0030728	0.00044206
% %     0.0024447	0.00020942	1.8759e-05
% %     0.00031798	1.3653e-05	6.2886e-07
% % ];
% % 
% % errorspost = [0.001552	0.00019875	2.4251e-05
% %     0.00010981	6.6628e-06	4.8759e-07
% %     7.9986e-06	2.3147e-07	1.6103e-08
% % ];

% clc
% clear all

% % When r=0.5
% % clc
% % errors=[0.10245 0.017186  0.0024447    %0.00031797
% %         0.031673  0.0030728  0.00020943  %1.3653e-05 
% %         0.0089015  0.00044206 1.8759e-05 ]; %6.2085e-07];
% %     
% %     
% % errorspost=[0.020924 0.001552   0.00010981   %7.9986e-06
% %             0.0046779  0.00019875  6.6628e-06  % 2.3147e-07  
% %             0.0011273  2.4251e-05  4.8759e-07 ];  %1.6103e-08];

%  % When r=0.41        
clc
errors=[ 0.29095      0.099231     0.027902  0.0072293
        0.10384  0.018087  0.0025258  0.00032583
        0.032316  0.0031196  0.00021821  1.4261e-05 
         0.008926  0.00046627  1.931e-05 6.3403e-07];
    
    
errorspost=[0.095842        0.022746     0.0054707  0.0013434
    0.020601  0.0015367  0.00011026  7.9099e-06
            0.0046738  0.00019954  6.785e-06  2.3502e-07
            0.0011272 2.4671e-05  4.9802e-07 8.1868e-09];


hs = 0.5.^[2 3 4 5];
lerrors = log10(errors);
lerrorspost = log10(errorspost);
lhs = log10(hs);


figure
semilogy(lhs (1,:),errors(1,:),'-ko','LineWidth',3,'MarkerSize',10)
hold on
semilogy(lhs (1,:),errors(2,:),'-ks','LineWidth',3,'MarkerSize',10)
hold on
semilogy(lhs (1,:),errors(3,:),'-k<','LineWidth',3,'MarkerSize',10)
hold on
semilogy(lhs (1,:),errors(4,:),'-k>','LineWidth',3,'MarkerSize',10)
hold on

semilogy(lhs (1,:),errorspost(1,:),'--ko','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
hold on
semilogy(lhs (1,:),errorspost(2,:),'--ks','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
hold on
semilogy(lhs (1,:),errorspost(3,:),'--k<','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
hold on
semilogy(lhs (1,:),errorspost(4,:),'--k>','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10)
hold on

[h]=legend('P1', 'P2','P3','P4','P1 postp.','P2 postp.','P3 postp.','P4 postp.','Location','southeast');
%h=findobj(ax,'type','axis');
%set(h(1), 'visible', 'off')    % Hides the legend's axes (legend border and background)
%set(h,'PlotBoxAspectRatioMode','manual');
%set(h,'PlotBoxAspectRatio',[1 0.8 1]);
set(h,'FontSize',20)
%set(icons,'FontSize',30)
%set(str,'FontSize',30)
xlabel('log_{10}(Mesh Size)','FontSize',20,'fontname','times','fontweight','bold'), ylabel('L2 error','FontSize',20,'fontname','times','fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold')
set(gca,'fontname','times')
%title('HDG','FontSize',15)
%legend boxoff 


p1slope1=(log10(errors(1,1))-log10(errors(1,2)))/(log10(0.5^2)-log10(0.5^3))
p1slope2=(log10(errors(1,2))-log10(errors(1,3)))/(log10(0.5^3)-log10(0.5^4))
p1slope3=(log10(errors(1,3))-log10(errors(1,4)))/(log10(0.5^4)-log10(0.5^5))

p2slope1=(log10(errors(2,1))-log10(errors(2,2)))/(log10(0.5^2)-log10(0.5^3))
p2slope2=(log10(errors(2,2))-log10(errors(2,3)))/(log10(0.5^3)-log10(0.5^4))
p2slope3=(log10(errors(2,3))-log10(errors(2,4)))/(log10(0.5^4)-log10(0.5^5))


p3slope1=(log10(errors(3,1))-log10(errors(3,2)))/(log10(0.5^2)-log10(0.5^3))
p3slope2=(log10(errors(3,2))-log10(errors(3,3)))/(log10(0.5^3)-log10(0.5^4))
p3slope3=(log10(errors(3,3))-log10(errors(3,4)))/(log10(0.5^4)-log10(0.5^5))

p4slope1=(log10(errors(4,1))-log10(errors(4,2)))/(log10(0.5^2)-log10(0.5^3))
p4slope2=(log10(errors(4,2))-log10(errors(4,3)))/(log10(0.5^3)-log10(0.5^4))
p4slope3=(log10(errors(4,3))-log10(errors(4,4)))/(log10(0.5^4)-log10(0.5^5))

p1postslope1=(log10(errorspost(1,1))-log10(errorspost(1,2)))/(log10(0.5^2)-log10(0.5^3))
p1postslope2=(log10(errorspost(1,2))-log10(errorspost(1,3)))/(log10(0.5^3)-log10(0.5^4))
p1postslope3=(log10(errorspost(1,3))-log10(errorspost(1,4)))/(log10(0.5^4)-log10(0.5^5))

p2postslope1=(log10(errorspost(2,1))-log10(errorspost(2,2)))/(log10(0.5^2)-log10(0.5^3))
p2postslope2=(log10(errorspost(2,2))-log10(errorspost(2,3)))/(log10(0.5^3)-log10(0.5^4))
p2postslope3=(log10(errorspost(2,3))-log10(errorspost(2,4)))/(log10(0.5^4)-log10(0.5^5))

p3postslope1=(log10(errorspost(3,1))-log10(errorspost(3,2)))/(log10(0.5^2)-log10(0.5^3))
p3postslope2=(log10(errorspost(3,2))-log10(errorspost(3,3)))/(log10(0.5^3)-log10(0.5^4))
p3postslope3=(log10(errorspost(3,3))-log10(errorspost(3,4)))/(log10(0.5^4)-log10(0.5^5))

p4postslope1=(log10(errorspost(4,1))-log10(errorspost(4,2)))/(log10(0.5^2)-log10(0.5^3))
p4postslope2=(log10(errorspost(4,2))-log10(errorspost(4,3)))/(log10(0.5^3)-log10(0.5^4))
p4postslope3=(log10(errorspost(4,3))-log10(errorspost(4,4)))/(log10(0.5^4)-log10(0.5^5))







% % hs = 0.5.^[3 4 5];
% % lhs = log10(hs);
% % figure (20)
% % 
% % semilogy(lhs (1,:),errors(1,:),'-ro','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errors(2,:),'-b>','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errors(3,:),'-g*','LineWidth',1)
% % 
% % semilogy(lhs (1,:),errorspost(1,:),'--ro','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errorspost(2,:),'--b>','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errorspost(3,:),'--g*','LineWidth',1)
% % 
% % legend('P2','P3','P4','P2 postprocesed','P3 postprocesed','P4 postprocesed')
% % xlabel('log_{10}(h)(Mesh Size)'), ylabel('L2 error')
% % title('Convergence X-HDG Dirichlet')




% slopesFirstSegment = (lerrorspost(end,:)-lerrorspost(end-1,:))/(lhs(end)-lhs(end-1));
% slopesSecondSegment = (lerrorspost(end-1,:)-lerrorspost(end-2,:))/(lhs(end-1)-lhs(end-2));
% slopesThirdSegment=(lerrors(end-2,:)-lerrors(end-3,:))/(lhs(end-2)-lhs(end-3))
% slopesFourthSegment=(lerrors(end-3,:)-lerrors(end-4,:))/(lhs(end-3)-lhs(end-4))

% errors2 = [9.7327e-16 0.0017094 0.0035983
%            1.8072e-15 0.0039668 0.0051079
%            2.3031e-15 0.0010137 0.0010137
%            2.2345e-15 0.0017646 0.0019422];
%        
% lerrors2 = log10(errors2);
% figure
% plot(lhs,lerrors2,'->')
% legend('P1','P2','P3')
% xlabel('log_{10}(h)'), ylabel('log_{10}(L2 error)')
% title('Convergence with New points')  %Pphy %This one actually makes no sense, the solution with new points is wrong!!!!!!!
% 
% 
% slopesFirstTwoPoints = (lerrors2(end,:)-lerrors2(end-1,:))/(lhs(end)-lhs(end-1))
% slopesLastTwoPoints = (lerrors2(end-1,:)-lerrors2(end-2,:))/(lhs(end-1)-lhs(end-2))
% slopesLastTwoPoints=(lerrors2(end-2,:)-lerrors2(end-3,:))/(lhs(end-2)-lhs(end-3))
% 
% errors3 = [9.7327e-16 1.2788e-09 7.3034e-06
%            1.8072e-15 5.6396e-06 0.00073849
%            2.3031e-15 6.6424e-08 5.0794e-05
%            2.2345e-15 5.529e-09  7.6356e-06
%            3.491e-15  2.6711e-10 9.5482e-07 ];
% lerrors3 = log10(errors3);
% figure
% plot(lhs,lerrors3,'-<')
% legend('P1','P2','P3')
% xlabel('log_{10}(h)'), ylabel('log_{10}(L2 error)')
% title('Convergence with New points (PtsInt)')  %PtsInt
% 
% 
% slopesFirstSegment = (lerrors(end,:)-lerrors(end-1,:))/(lhs(end)-lhs(end-1))
% slopesSecondSegment = (lerrors(end-1,:)-lerrors(end-2,:))/(lhs(end-1)-lhs(end-2))
% slopesThirdSegment=(lerrors(end-2,:)-lerrors(end-3,:))/(lhs(end-2)-lhs(end-3))
% slopesFourthSegment=(lerrors(end-3,:)-lerrors(end-4,:))/(lhs(end-3)-lhs(end-4))
% 
% 
%  SlopeTriangles (lhs,lerrorspost,slopesFirstSegment);



















