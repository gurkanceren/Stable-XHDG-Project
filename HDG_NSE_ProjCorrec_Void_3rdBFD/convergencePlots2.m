
function convergencePlots2

clc

           
%%%%% with simple analytical solution single material pure HDG
% % 
% % errors=[0.13973  0.018168  0.0022938 
% %          0.0091789   0.00058227  3.6517e-05 
% %           0.0003011  9.238e-06   2.8608e-07  ];
% % 
% % 
% % 
% % errorspost=[ 0.03599  0.0030959  0.00026131 
% %               0.0023507  0.00010205  4.4373e-06 
% %                 7.6007e-05   1.6313e-06   3.5495e-08];
            
                       
%%%%% with simple analytical solution single material pure HDG tau=1 single
%%%%% face

% % errors=[0.44473  0.05301  0.006427  
% %          0.024735    0.0014538  8.894e-05 
% %            0.00039183   1.2075e-05  3.7755e-07  ];
% % 
% % 
% % 
% % errorspost=[0.032833   0.0030004  0.00025832 
% %              0.0021859    9.8468e-05   4.3582e-06 
% %               7.1539e-05   1.5808e-06  3.4931e-08  ];

         
         
hs = 0.5.^[2 3 4];
lerrors = log10(errors);
lerrorspost = log10(errorspost);
lhs = log10(hs);
figure (22)


figure
semilogy(lhs (1,:),errors(1,:),'-bo','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost(1,:),'--bo','LineWidth',3,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',12)
hold on


semilogy(lhs (1,:),errors(2,:),'-rs','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost(2,:),'--rs','LineWidth',3,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12)
hold on


semilogy(lhs (1,:),errors(3,:),'-g<','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost(3,:),'--g<','LineWidth',3,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',12)
hold on


% semilogy(lhs (1,:),errors(4,:),'-md','LineWidth',3,'MarkerSize',12)
% hold on
% semilogy(lhs (1,:),errorspost(4,:),'--md','LineWidth',3,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',12)
% hold on


[h]=legend('P2','P2 postp.','P3','P3 postp.','P4','P4 postp.','Location','southeast');
%h=findobj(ax,'type','axis');
%set(h(1), 'visible', 'off')    % Hides the legend's axes (legend border and background)
set(h,'PlotBoxAspectRatioMode','manual');
set(h,'PlotBoxAspectRatio',[1 1.1 1]);
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
%p1slope3=(log10(errors(1,3))-log10(errors(1,4)))/(log10(0.5^3)-log10(0.5^4))

p2slope1=(log10(errors(2,1))-log10(errors(2,2)))/(log10(0.5^1)-log10(0.5^2))
p2slope2=(log10(errors(2,2))-log10(errors(2,3)))/(log10(0.5^2)-log10(0.5^3))
%p2slope3=(log10(errors(2,3))-log10(errors(2,4)))/(log10(0.5^3)-log10(0.5^4))

p3slope1=(log10(errors(3,1))-log10(errors(3,2)))/(log10(0.5^1)-log10(0.5^2))
p3slope2=(log10(errors(3,2))-log10(errors(3,3)))/(log10(0.5^2)-log10(0.5^3))
%p3slope3=(log10(errors(3,3))-log10(errors(3,4)))/(log10(0.5^3)-log10(0.5^4))

% p4slope1=(log10(errors(4,1))-log10(errors(4,2)))/(log10(0.5^1)-log10(0.5^2))
% p4slope2=(log10(errors(4,2))-log10(errors(4,3)))/(log10(0.5^2)-log10(0.5^3))
% p4slope3=(log10(errors(4,3))-log10(errors(4,4)))/(log10(0.5^3)-log10(0.5^4))


p1postslope1=(log10(errorspost(1,1))-log10(errorspost(1,2)))/(log10(0.5^1)-log10(0.5^2))
p1postslope2=(log10(errorspost(1,2))-log10(errorspost(1,3)))/(log10(0.5^2)-log10(0.5^3))
%p1postslope3=(log10(errorspost(1,3))-log10(errorspost(1,4)))/(log10(0.5^3)-log10(0.5^4))

p2postslope1=(log10(errorspost(2,1))-log10(errorspost(2,2)))/(log10(0.5^1)-log10(0.5^2))
p2postslope2=(log10(errorspost(2,2))-log10(errorspost(2,3)))/(log10(0.5^2)-log10(0.5^3))
%p2postslope3=(log10(errorspost(2,3))-log10(errorspost(2,4)))/(log10(0.5^3)-log10(0.5^4))

p3postslope1=(log10(errorspost(3,1))-log10(errorspost(3,2)))/(log10(0.5^1)-log10(0.5^2))
p3postslope2=(log10(errorspost(3,2))-log10(errorspost(3,3)))/(log10(0.5^2)-log10(0.5^3))
%p3postslope3=(log10(errorspost(3,3))-log10(errorspost(3,4)))/(log10(0.5^3)-log10(0.5^4))

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
 