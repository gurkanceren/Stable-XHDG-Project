

function convergencePlots

errors = [9.7327e-16	1.6176e-15	9.9463e-09
    1.8072e-15	1.8795e-15	1.6944e-06
    2.3031e-15	6.8739e-15	2.2959e-07
    2.2345e-15  1.1722e-14  1.3541e-08
    3.491e-15   3.7417e-14  5.4163e-10];

hs = 0.5.^[1 2 3 4 5];
lerrors = log10(errors);
lhs = log10(hs);
figure
plot(lhs,lerrors,'o-')
legend('P1','P2','P3')
xlabel('log_{10}(h)'), ylabel('log_{10}(L2 error)')
title('Convergence with Esther points')


slopesFirstSegment = (lerrors(end,:)-lerrors(end-1,:))/(lhs(end)-lhs(end-1))
slopesSecondSegment = (lerrors(end-1,:)-lerrors(end-2,:))/(lhs(end-1)-lhs(end-2))
slopesThirdSegment=(lerrors(end-2,:)-lerrors(end-3,:))/(lhs(end-2)-lhs(end-3))
slopesFourthSegment=(lerrors(end-3,:)-lerrors(end-4,:))/(lhs(end-3)-lhs(end-4))

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
errors3 = [9.7327e-16 1.2788e-09 7.3034e-06
           1.8072e-15 5.6396e-06 0.00073849
           2.3031e-15 6.6424e-08 5.0794e-05
           2.2345e-15 5.529e-09  7.6356e-06
           3.491e-15  2.6711e-10 9.5482e-07 ];
lerrors3 = log10(errors3);
figure
plot(lhs,lerrors3,'-<')
legend('P1','P2','P3')
xlabel('log_{10}(h)'), ylabel('log_{10}(L2 error)')
title('Convergence with New points (PtsInt)')  %PtsInt


slopesFirstSegment = (lerrors(end,:)-lerrors(end-1,:))/(lhs(end)-lhs(end-1))
slopesSecondSegment = (lerrors(end-1,:)-lerrors(end-2,:))/(lhs(end-1)-lhs(end-2))
slopesThirdSegment=(lerrors(end-2,:)-lerrors(end-3,:))/(lhs(end-2)-lhs(end-3))
slopesFourthSegment=(lerrors(end-3,:)-lerrors(end-4,:))/(lhs(end-3)-lhs(end-4))


SlopeTriangles (lhs,lerrors3,slopesFirstSegment);



















