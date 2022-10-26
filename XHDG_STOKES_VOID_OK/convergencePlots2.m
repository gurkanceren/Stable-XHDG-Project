
function convergencePlots2

%% Void Stokes Neumann convergence 
%errors=[0.038686  0.005064  0.00064361  8.0786e-05 
%          0.0024138  0.00015408  9.8132e-06  6.1829e-07 
%           7.3307e-05  2.299e-06  7.3145e-08  2.3222e-09 ] ;
% 
% 
% 
% errorspost=[0.0034473  0.00023639  1.5549e-05  9.8886e-07 
%              0.00015748   5.1044e-06  1.65e-07  5.2383e-09 
%                2.9548e-06  4.6304e-08  7.4176e-10  7.6507e-11 ];
  
%% Void Stokes Neumann convergence without mesh modification
% errors=[ 0.0386    0.0050  6.3783e-04   8.0110e-05
%          0.0024  1.5196e-04   9.6020e-06   6.0411e-07
%          7.1767e-05   2.2502e-06  7.0841e-08  2.7815e-06 ] ;
% 
% 
% 
% errorspost=[0.0034   2.3556e-04   1.5365e-05   9.8009e-07 
%        1.5266e-04    5.0096e-06   1.6123e-07   5.1202e-09 
%         2.9472e-06  5.4405e-08   1.2628e-09 2.2868e-07];

    
%% Void Stokes Neumann convergence with mesh modification
errors=[ 0.038558    0.0050339   0.00063783    8.0236e-05 
          0.0023841  0.00015196  9.602e-06  6.0601e-07 
            7.1767e-05   2.2502e-06  7.0841e-08  2.2368e-09 ] ;



errorspost=[0.0034134   0.00023556    1.5365e-05   9.8792e-07
            0.00015266  5.0096e-06  1.6123e-07  5.1307e-09 
             2.9472e-06  5.4405e-08  1.2628e-09   1.1949e-11 ];
         
         
hs = 0.5.^[1 2 3 4];
lerrors = log10(errors);
lerrorspost = log10(errorspost);
lhs = log10(hs);
figure
x0=10;
y0=10;
width=500;
height=400;
set(gcf,'units','points','position',[x0,y0,width,height])
ylim([1e-5 10^-4])
xlim([10^0.2  10^1.2])

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

[h]=legend('P2','P2 postp.','P3','P3 postp.','P4','P4 postp.','Location','southeast');
%[h]=legend('P2','P2 postp.','P3','P3 postp.','P4','P4 postp.','Location','southeast');
%h=findobj(ax,'type','axis');
%set(h(1), 'visible', 'off')    % Hides the legend's axes (legend border and background)
%set(h,'PlotBoxAspectRatioMode','manual');
%set(h,'PlotBoxAspectRatio',[1 0.8 1]);
%set(h,'FontSize',20)
%set(icons,'FontSize',30)
%set(str,'FontSize',30)
xlabel('log_{10}(Mesh Size)','FontSize',20,'fontname','times','fontweight','bold'), ylabel('L2 error','FontSize',20,'fontname','times','fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')
set(gca,'fontname','times')
ylim([1e-11 1e-1])
xlim([-1.3  -0.3])
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
 