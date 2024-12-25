

function convergencePlots_test

% % errors = [0.017186	0.0030728	0.00044206
% %     0.0024447	0.00020942	1.8759e-05
% %     0.00031798	1.3653e-05	6.2886e-07
% % ];
% % 
% % errorspost = [0.001552	0.00019875	2.4251e-05
% %     0.00010981	6.6628e-06	4.8759e-07
% %     7.9986e-06	2.3147e-07	1.6103e-08
% % ];

clc
clear all
errors=[0.10245 0.017186  0.0024447
        0.0089015  0.00044206 1.8759e-05];
    
    
errorspost=[0.020924 0.001552   0.00010981
            0.0011273  2.4251e-05  4.8759e-07];


hs = 0.5.^[2 3 4];
lerrors = log10(errors);
lerrorspost = log10(errorspost);
lhs = log10(hs);


figure
semilogy(lhs (1,:),errors(1,:),'-ko','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errors(2,:),'-k>','LineWidth',3,'MarkerSize',12)
hold on

semilogy(lhs (1,:),errorspost(1,:),'--ko','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost(2,:),'--k>','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12)
hold on

[h]=legend('P2','P4','P2 postp.','P4 postp.','Location','southeast');
%leg=findobj(ax,'type','text');
set(h,'PlotBoxAspectRatioMode','manual');
set(h,'PlotBoxAspectRatio',[1 0.6 1]);
set(h,'FontSize',30)
%set(icons,'FontSize',30)
%set(str,'FontSize',30)
xlabel('log_{10}(Mesh Size)','FontSize',30), ylabel('L2 error','FontSize',30)
%title('HDG','FontSize',15)


p2slope1=(log10(errors(1,1))-log10(errors(1,2)))/(log10(0.5^2)-log10(0.5^3))

p2slope2=(log10(errors(1,2))-log10(errors(1,3)))/(log10(0.5^3)-log10(0.5^4))


p4slope1=(log10(errors(2,1))-log10(errors(2,2)))/(log10(0.5^2)-log10(0.5^3))

p4slope2=(log10(errors(2,2))-log10(errors(2,3)))/(log10(0.5^3)-log10(0.5^4))


p2postslope1=(log10(errorspost(1,1))-log10(errorspost(1,2)))/(log10(0.5^2)-log10(0.5^3))

p2postslope2=(log10(errorspost(1,2))-log10(errorspost(1,3)))/(log10(0.5^3)-log10(0.5^4))


p4postslope1=(log10(errorspost(2,1))-log10(errorspost(2,2)))/(log10(0.5^2)-log10(0.5^3))

p4postslope2=(log10(errorspost(2,2))-log10(errorspost(2,3)))/(log10(0.5^3)-log10(0.5^4))










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



















