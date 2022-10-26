
 

figure
semilogy(lhs (1,:),errors1(1,:),'-bo','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost1(1,:),'--bo','LineWidth',3,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',12)
hold on


semilogy(lhs (1,:),errors1(2,:),'-rs','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost1(2,:),'--rs','LineWidth',3,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12)
hold on


semilogy(lhs (1,:),errors1(3,:),'-g<','LineWidth',3,'MarkerSize',12)
hold on
semilogy(lhs (1,:),errorspost1(3,:),'--g<','LineWidth',3,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',12)
hold on



[h]=legend('P1','P1 postp.','P2','P2 postp.','P3','P3 postp.','Location','southeast');
%h=findobj(ax,'type','axis');
%set(h(1), 'visible', 'off')    % Hides the legend's axes (legend border and background)
set(h,'PlotBoxAspectRatioMode','manual');
set(h,'PlotBoxAspectRatio',[1 1.1 1]);
set(h,'FontSize',15)
%set(icons,'FontSize',30)
%set(str,'FontSize',30)
xlabel('log_{10}(Mesh Size)','FontSize',20,'fontname','times','fontweight','bold'), ylabel('L2 error','FontSize',20,'fontname','times','fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')
set(gca,'fontname','times')
%title('HDG','FontSize',15)
%legend boxoff 

p2slope1=(log10(errors1(1,1))-log10(errors1(1,2)))/(log10(0.5^1)-log10(0.5^2))
p2slope2=(log10(errors1(1,2))-log10(errors1(1,3)))/(log10(0.5^2)-log10(0.5^3))
%p2slope3=(log10(errors1(1,3))-log10(errors1(1,4)))/(log10(0.5^3)-log10(0.5^4))

p3slope1=(log10(errors1(2,1))-log10(errors1(2,2)))/(log10(0.5^1)-log10(0.5^2))
p3slope2=(log10(errors1(2,2))-log10(errors1(2,3)))/(log10(0.5^2)-log10(0.5^3))
%p3slope3=(log10(errors1(2,3))-log10(errors1(2,4)))/(log10(0.5^3)-log10(0.5^4))

p4slope1=(log10(errors1(3,1))-log10(errors1(3,2)))/(log10(0.5^1)-log10(0.5^2))
p4slope2=(log10(errors1(3,2))-log10(errors1(3,3)))/(log10(0.5^2)-log10(0.5^3))
%p4slope3=(log10(errors1(3,3))-log10(errors1(3,4)))/(log10(0.5^3)-log10(0.5^4))

% % p5slope1=(log10(errors(4,1))-log10(errors(4,2)))/(log10(0.5^1)-log10(0.5^2))
% % p5slope2=(log10(errors(4,2))-log10(errors(4,3)))/(log10(0.5^2)-log10(0.5^3))
% % p5slope3=(log10(errors(4,3))-log10(errors(4,4)))/(log10(0.5^3)-log10(0.5^4))

p2postslope1=(log10(errorspost1(1,1))-log10(errorspost1(1,2)))/(log10(0.5^1)-log10(0.5^2))
p2postslope2=(log10(errorspost1(1,2))-log10(errorspost1(1,3)))/(log10(0.5^2)-log10(0.5^3))
%p2postslope3=(log10(errorspost1(1,3))-log10(errorspost1(1,4)))/(log10(0.5^3)-log10(0.5^4))

p3postslope1=(log10(errorspost1(2,1))-log10(errorspost1(2,2)))/(log10(0.5^1)-log10(0.5^2))
p3postslope2=(log10(errorspost1(2,2))-log10(errorspost1(2,3)))/(log10(0.5^2)-log10(0.5^3))
%p3postslope3=(log10(errorspost1(2,3))-log10(errorspost1(2,4)))/(log10(0.5^3)-log10(0.5^4))

p4postslope1=(log10(errorspost1(3,1))-log10(errorspost1(3,2)))/(log10(0.5^1)-log10(0.5^2))
p4postslope2=(log10(errorspost1(3,2))-log10(errorspost1(3,3)))/(log10(0.5^2)-log10(0.5^3))
%p4postslope3=(log10(errorspost1(3,3))-log10(errorspost1(3,4)))/(log10(0.5^3)-log10(0.5^4))

% % p5postslope1=(log10(errorspost(4,1))-log10(errorspost(4,2)))/(log10(0.5^1)-log10(0.5^2))
% % p5postslope2=(log10(errorspost(4,2))-log10(errorspost(4,3)))/(log10(0.5^2)-log10(0.5^3))
% % p5postslope3=(log10(errorspost(4,3))-log10(errorspost(4,4)))/(log10(0.5^3)-log10(0.5^4))




