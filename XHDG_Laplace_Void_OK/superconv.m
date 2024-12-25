
function superconv

errors = [0.00020943	1.87e-05  1.73e-06];

errorspost = [6.66e-06	4.87e-07	2.80e-07];


hs = 0.5.^[4];
lerrors = log10(errors);
lerrorspost = log10(errorspost);
lhs = log10(hs);
figure
plot(lhs,lerrors,'ro-','LineWidth',1)
hold on
plot(lhs,lerrorspost,'b<--','LineWidth',1.2)

legend('P3','P4','P5','P3 postprocesed','P4 postprocesed','P5 postprocesed')
xlabel('log_{10}(h)(Mesh Size)'), ylabel('log_{10}(L2 error)')
title('Convergence HDG')




% slopesFirstSegment = (lerrorspost(end,:)-lerrorspost(end-1,:))/(lhs(end)-lhs(end-1));
% slopesSecondSegment = (lerrorspost(end-1,:)-lerrorspost(end-2,:))/(lhs(end-1)-lhs(end-2));
% 
% SlopeTriangles (lhs,lerrors,slopesFirstSegment);