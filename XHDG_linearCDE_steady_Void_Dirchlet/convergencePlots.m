

function convergencePlots

errors = [0.017186	0.0030727	0.00044206
    0.0024447	0.00020942	1.8759e-05
    0.00031798	1.3653e-05	6.2886e-07
];

errorspost = [0.001552	0.00019875	2.4251e-05
    0.00010981	6.6628e-06	4.8759e-07
    7.9986e-06	2.3147e-07	1.6103e-08
];

% % hs = 0.5.^[3 4 5];
% % lerrors = log10(errors);
% % lerrorspost = log10(errorspost);
% % lhs = log10(hs);
% % 
% % figure
% % hold on
% % plot(lhs (1,:),lerrorspost(1,:),'-ko','LineWidth',1.1)
% % hold on
% % plot(lhs (1,:),lerrorspost(2,:),'-k>','LineWidth',1.1)
% % hold on
% % plot(lhs (1,:),lerrorspost(3,:),'-ks','LineWidth',1.1)
% % 
% % 
% % legend('P2 postprocessed','P3 postprocessed','P4 postprocessed')
% % 
% % 
% % xlabel('log_{10}(h)(Mesh Size)'), ylabel('log_{10}(L2 error)')
% % title('X-HDG Postprocessed Convergence')


        
        
        
hs = 0.5.^[2 3 4];
lhs = log10(hs);
figure (20)

semilogy(lhs (1,:),errors(1,:),'-ko','LineWidth',1)
hold on
semilogy(lhs (1,:),errors(2,:),'-k>','LineWidth',1)
hold on
semilogy(lhs (1,:),errors(3,:),'-ks','LineWidth',1)

semilogy(lhs (1,:),errorspost(1,:),'--ko','LineWidth',1)
hold on
semilogy(lhs (1,:),errorspost(2,:),'--k>','LineWidth',1)
hold on
semilogy(lhs (1,:),errorspost(3,:),'--ks','LineWidth',1)

legend('P2','P3','P4','P2 postprocesed','P3 postprocesed','P4 postprocesed')
xlabel('log_{10}(h)(Mesh Size)'), ylabel('L2 error')
title('Convergence X-HDG')


%  slopesFirstSegment = (lerrorspost(end,:)-lerrorspost(end-1,:))/(lhs(end)-lhs(end-1));
%  slopesFirstSegment = ceil(slopesFirstSegment*100)/100;
%  SlopeTriangles (lhs,lerrorspost,slopesFirstSegment);



















