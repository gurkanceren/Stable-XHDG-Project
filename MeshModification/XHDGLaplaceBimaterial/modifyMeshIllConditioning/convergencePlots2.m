
function convergencePlots2

errors = [0.017186	0.0030727	0.00044206
    0.0024447	0.00020942	1.8759e-05
    0.00031798	1.3653e-05	6.2886e-07
];

errorspost = [0.003155 0.0025  0.0025
              0.00041  0.000396  0.000396
              7.96e-05  7.92e-05  7.92e-05];


hs = 0.5.^[3 4 5];
lerrors = log10(errors);
lerrorspost = log10(errorspost);
lhs = log10(hs);
figure
% plot(lhs,lerrors,'ro-','LineWidth',1)
% hold on
plot(lhs,lerrorspost,'b<--','LineWidth',1.2)

legend('P2','P3','P4','P2 postprocesed','P3 postprocesed','P4 postprocesed')
xlabel('log_{10}(h)(Mesh Size)'), ylabel('log_{10}(L2 error)')
title('Convergence HDG')



% % 
% % hs = 0.5.^[3 4 5];
% % lhs = log10(hs);
% % figure (20)
% % 
% % semilogy(lhs (1,:),errors(1,:),'-ro','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errors(2,:),'-b>','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errors(3,:),'-g<','LineWidth',1)
% % 
% % semilogy(lhs (1,:),errorspost(1,:),'-ro','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errorspost(2,:),'-b>','LineWidth',1)
% % hold on
% % semilogy(lhs (1,:),errorspost(3,:),'-g<','LineWidth',1)
% % 
% % legend('P2','P3','P4','P2 postprocesed','P3 postprocesed','P4 postprocesed')
% % xlabel('log_{10}(h)(Mesh Size)'), ylabel('L2 error')
% % title('Convergence HDG')




slopesFirstSegment = (lerrorspost(end,:)-lerrorspost(end-1,:))/(lhs(end)-lhs(end-1));
slopesSecondSegment = (lerrorspost(end-1,:)-lerrorspost(end-2,:))/(lhs(end-1)-lhs(end-2));
 SlopeTriangles (lhs,lerrorspost,slopesFirstSegment);
 
 