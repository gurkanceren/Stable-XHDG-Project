

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

% clc
% errors=[0.12703  0.033203  0.0083383
%         0.0050554  0.00068354  8.4627e-05 
%          2.2896e-14   1.2462e-12  2.6792e-14 
%           8.7854e-14   9.0695e-11  1.8113e-13 ];
%     
%     
% errorspost=[ 0.08286  0.01998   0.0048694
%              0.0022305  0.00019082  1.6696e-05 
%              9.1186e-15  2.4952e-13 1.7243e-14 
%               2.2385e-14  1.891e-11 4.767e-14 ];
% clc

%%%% Diffrent analytical soln x=0.2 poly order 6 !! 
% % errors=[5.027  1.7356  0.47061 0.11963
% %          1.1234   0.19214  0.025706 0.0032433
% %           0.12306   0.01115  0.00077427 4.9792e-05
% %            0.0070044   0.00031267 1.6054e-05  1.1254e-05];
% %     
% %     
% % errorspost=[2.0605  0.66165 0.15864 0.037912
% %             0.3708  0.050392  0.0043702 0.00036934
% %             0.066323   0.0032563  0.0001441  1.2878e-05
% %              0.005374   0.00010426  1.1488e-05 1.1248e-05];

clc        
%%%% Diffrent analytical soln          
% % errors=[13.9833  5.1477 1.4696 0.37861 
% %          4.1408  0.79377 0.10867 0.013874 
% %           0.69202 0.070528  0.0047589  0.00030091
% %              0.060494  0.0032747  0.00011537  1.0086e-05];
% %     
% %     
% % errorspost=[4.5401  1.5859  0.38935  0.092293 
% %              1.0033  0.16373  0.015054 0.0012763 
% %               0.22751  0.014776  0.00067043 3.0804e-05
% %                 0.033176  0.00075169  1.8727e-05  9.4254e-06];
            
 
%%% x=0.1875 poly order 6 with high diffreence in viscosity !!
% % errors=[27.9061  10.2825   2.9361  0.75646 
% %         8.4315   1.5867  0.21713  0.030074 
% %           1.4372  0.14105  0.0095089   0.00060099  
% %            0.12672  0.0065501  0.00022977  7.9126e-06];
% %     
% %     
% % errorspost=[9.1131  3.1669  0.77746   0.18428
% %             2.0231   0.32699  0.030074 0.0025481
% %              0.45863  0.029531  0.0013393  5.8622e-05 
% %               0.066654  0.0015026  3.2676e-05 3.3384e-06 ];
          
%%% x=0.1875 poly order 6 with high diffreence in viscosity !!
% % errors=[ 
% %         8.4315   1.5867  0.21713  
% %           1.4372  0.14105  0.0095089   
% %            0.12672  0.0065501  0.00022974 ];
% %     
% %     
% % errorspost=[
% %             2.0231   0.32699  0.030074
% %              0.45863  0.029532  0.0013393   
% %               0.066654  0.0015027  3.2485e-05 ];


% % %%% x=0.1875 poly order 6 with slight difference in vscosity
% % 
% % errors=[40.3146  14.2508 4.028 1.0354 
% %          11.8309  2.1699  0.29585  0.037761
% %           2.3001  0.19525  0.012977 0.00082055
% %            0.2491  0.009704 0.00032441 1.0295e-05 ];
% %     
% %     
% % errorspost=[15.3485  4.6059 1.1141  0.26503
% %             3.3204  0.46114  0.041758  0.0035089
% %              0.68902  0.040426  0.0018247 7.9745e-05 
% %               0.095479  0.0020882 4.4578e-05  9.6285e-07];
% %           
%%% x=0.1875 poly order 6 with slight difference in vscosity

% % errors=[
% %          11.8309  2.1699  0.29585  0.037761
% %           2.3001  0.19525  0.012977 0.00082055
% %            0.2491  0.009704 0.00032441 1.0295e-05 ];
% %     
% %     
% % errorspost=[
% %             3.3204  0.46114  0.041758  0.0035089
% %              0.68902  0.040426  0.0018247 7.9745e-05 
% %               0.095479  0.0020882 4.4578e-05  9.6285e-07];
          

% % errors=[0.17438  0.027527  0.0036272  0.00045776
% %          0.0232  0.0017092  0.0001135   7.182e-06
% %           0.0015927  5.2962e-05   1.7391e-06  5.5443e-08 ];
% %     
% %     
% % errorspost=[0.059139  0.0071835   0.0006191  5.226e-05
% %              0.01013  0.00046841 2.0387e-05  8.8807e-07  
% %               0.00076848   1.605e-05  3.3263e-07  7.2147e-09 ];


%%%% final version for linear interface bimaterial problem polynomial
%%%% solution mu1=1 mu2=2.5
          
% % errors=[2.684  0.94011  0.25348  0.064356
% %         0.69426  0.10598  0.013844  0.001747  
% %           0.10018  0.0068613  0.00044767  2.8325e-05 
% %             0.0072442  0.00022304  7.0147e-06   2.259e-07];
% %     
% %     
% % errorspost=[1.2316   0.36564  0.087  0.020829  
% %              0.23848  0.027503  0.0023623  0.00019925  
% %               0.040638  0.0017991    7.7978e-05  3.3961e-06 
% %                 0.0029085  6.3772e-05  1.3595e-06 2.9135e-08 ];
            
            
%%%% final version for curved interface bimaterial problem polynomial
%%%% solution mu1=2 mu2=1           


% % errors=[9.2186  2.6558  0.73702  0.19168  
% %          1.6841  0.35726  0.048844  0.0063078 
% %            0.26593  0.033893  0.0022477 0.00014487 
% %              0.02649   0.0018674  6.3e-05   4.405e-06];
% %     
% %     
% % errorspost=[7.1205   1.7106  0.41437  0.1016 
% %              0.78467   0.087932  0.0068681  0.00053648 
% %               0.10422  0.0055859   0.00022173   8.9483e-06
% %                 0.0082315   0.00020303    3.9602e-06  7.9955e-08 ];

            
%%%% final version for curved interface bimaterial problem complex
%%%% solution mu1=1 mu2=5  r=0.4         


% % errors=[0.87123  0.28066  0.086467  0.023107 
% %          0.31299  0.061476  0.010205   0.0013776 
% %            0.045884  0.016716  0.0014016  9.3085e-05 
% %              0.038106  0.0042785  0.00014758  0.0014143  ];  
% %     
% %     
% % errorspost=[0.29871 0.11152  0.024698  0.0058929 
% %              0.1222  0.014426  0.00083575  5.9597e-05 
% %               0.019616  0.0021228  9.0385e-05   3.0455e-06 
% %                  0.0047207  0.00050451  7.4932e-06  2.5465e-06 ];  
             
             
%%%% final version for curved interface bimaterial problem complex
%%%% solution mu1=1 mu2=5  r=0.41         


% % errors=[ m3p4 0.00014756  m4p4 7.4838e-06]; 
% %     
% %     
% % errorspost= [m3p4 7.6003e-06  m4p4 1.8063e-06  ];
             

%%%% final version for curved interface bimaterial problem complex
%%%% solution mu1=1 mu2=5  r=0.4125         


% % errors=[ 0.30573   0.061382   0.010205    0.0013776 
% %          0.04464   0.016716   0.0014013   9.2741e-05 
% %          0.037589  0.0042772  0.00014754  5.4401e-06 ]; 
% %     
% %     
% % errorspost= [  
% %           0.12114    0.014123    0.00083569   5.9327e-05    
% %           0.019246   0.0021366   9.0289e-05   2.9869e-06  
% %           0.0045992  0.00049361   7.6546e-06  1.3832e-07 ];

%%%% final version for curved interface bimaterial problem complex
%%%% solution mu1=1 mu2=5  r=0.425         

% % 
% % errors=[m4p4 5.4337e-06 m5p4 1.1164e-05  ]; 
% %     
% %     
% % errorspost= [m4p4 1.3547e-07  m5p4   1.4696e-07 ];

%%%% final version for curved interface bimaterial problem complex
%%%% solution mu1=1 mu2=5  r=0.43        


% % errors=[  0.062307   0.010204    0.0013776    0.00017581  
% %           0.016717   0.0014008   9.2728e-05   5.92e-06  
% %           0.0042757  0.00014752  5.4337e-06   1.9702e-07 ]; 
% %     
% %     
% % errorspost= [  0.017088    0.00083561  5.9783e-05   4.5105e-06  
% %                0.0021659   9.0271e-05  2.9934e-06   1.1323e-07 
% %                0.00048016  7.7556e-06  1.3568e-07   5.8723e-09 ];
           
           
%%%% final version for curved interface bimaterial problem complex
%%%% solution mu1=1 mu2=5  r=0.43        


% % errors=[  0.062307   0.010204    0.0013776      
% %           0.016717   0.0014008   9.2728e-05    
% %           0.0042757  0.00014752  5.4337e-06 ]; 
% %     
% %     
% % errorspost= [  0.017088    0.00083561  5.9783e-05   
% %                0.0021659   9.0271e-05  2.9934e-06  
% %                0.00048016  7.7556e-06  1.3568e-07];
             



%% final version for curved interface bimaterial problem 
%% solution mu1=1 mu2=100  r=0.5        

% 
% errors=[ 0.6225    0.1682    0.0430    0.0108  
%          0.0751    0.0098    0.0013    0.0002
%          0.00490102419560172	0.000310401938348694	1.96295981934980e-05	1.24125993996462e-06]; %%row 1 means P1 m1,m2,m3,m4
%     
%     
% errorspost= [ 0.0747    0.0117    0.0016    0.0002 
%            0.00758355286084166	0.000489628256672612	3.11046241387079e-05	1.96968155581470e-06
%            0.000467353169125251	1.60575226298701e-05	5.00578130723136e-07	1.68362456668330e-08];

%% laplace void problem         

% 
% errors=[ 0.10245 0.017186  0.0024447  0.00031797 
%          0.031673  3.0728e-3  2.0943e-4  1.3653e-5 
%            0.0089015  4.4208e-4  1.8759e-5  6.2118e-7]; %%row 1 means P1 m1,m2,m3,m4
%     
%     
% errorspost= [0.020924  0.001552 0.00010981   7.9986e-6 
%               0.0046779 1.9875e-4  6.6628e-6  2.3147e-7 
%                0.0011273  2.425e-5  4.8759e-7  2.4219e-8];

% %% stokes bimat prob
% errors=[ 0.9040    0.3892    0.1975    0.0536
%         0.2811  0.1398   0.0051  6.4508e-04
%          0.0217   0.0024  1.5592e-04   9.9888e-06 ] ;
% 
% 
% 
% errorspost=[ 0.2813    0.1392    0.0058    0.0008
%            0.0341   0.0035   2.3644e-04   1.5557e-05
%              0.0038  1.5763e-04  5.1704e-06   2.0426e-07 ];

%% stokes void prob
% errors=[ 0.15249 0.038558    0.0050339   0.00063783  
%           0.019541 0.0023841  0.00015196  9.602e-06  
%             0.0015841  7.1767e-05   2.2502e-06  7.0841e-08  ] ;
% 
% 
% 
% errorspost=[0.029315 0.0034134   0.00023556    1.5365e-05  
%             0.0034312  0.00015266  5.0096e-06  1.6123e-07  
%             0.00011415  2.9472e-06  5.4405e-08  1.2628e-09  ];

%% Laplace void problem


% errors = [ 0.10245  0.017186	0.0030727	0.00044206
%              0.031675  0.0024447	0.00020942	1.8759e-05
%                    0.0089008   0.00031798	1.3653e-05	6.2886e-07 ];


errors = [0.28882  0.10245   0.017186  0.0024447
0.17945    0.031675  0.0030727  0.00020942
0.073188   0.0089008   0.00044206  1.8759e-05 ] ;





  
% errorspost = [0.020921  0.001552	0.00019875	2.4251e-05
%    0.0046815  0.00010981	6.6628e-06	4.8759e-07
%     0.0011305  7.9986e-06	2.3147e-07	1.6103e-08];

    
errorspost = [ 0.1092     0.020921   0.001552    0.00010981
0.043253   0.0046815  0.00019879   6.6636e-06
0.013436   0.0011305   2.4258e-05   4.8682e-07 ];
             
  
hs = 0.5.^[0 1 2 3];
lerrors = log10(errors);
lerrorspost = log10(errorspost);
lhs = log10(hs);


figure
semilogy(lhs (1,:),errors(1,:),'-bo','LineWidth',2,'MarkerSize',6)
hold on
semilogy(lhs (1,:),errorspost(1,:),'--bo','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
hold on


semilogy(lhs (1,:),errors(2,:),'-rs','LineWidth',2,'MarkerSize',6)
hold on
semilogy(lhs (1,:),errorspost(2,:),'--rs','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6)
hold on


semilogy(lhs (1,:),errors(3,:),'-g<','LineWidth',2,'MarkerSize',6)
hold on
semilogy(lhs (1,:),errorspost(3,:),'--g<','LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
hold on


% semilogy(lhs (1,:),errors(4,:),'-md','LineWidth',3,'MarkerSize',12)
% hold on
% semilogy(lhs (1,:),errorspost(4,:),'--md','LineWidth',3,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',12)
% hold on

[h]=legend('P2','P2 postp.','P3','P3 postp.','P4','P4 postp.','Location','southeast');
%h=findobj(ax,'type','axis');
%set(h(1), 'visible', 'off')    % Hides the legend's axes (legend border and background)
%set(h,'PlotBoxAspectRatioMode','manual');
%set(h,'PlotBoxAspectRatio',[1.5 1.5 1.5]);
set(h,'FontSize',12)
%set(icons,'FontSize',30)
%set(str,'FontSize',30)
xlabel('log_{10}(Mesh Size)','FontSize',12,'fontname','times','fontweight','bold'), ylabel('L2 error','FontSize',12,'fontname','times','fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold')
set(gca,'fontname','times')
%title('HDG','FontSize',15)
%legend boxoff 


p2slope1=(log10(errors(1,1))-log10(errors(1,2)))/(log10(0.5^1)-log10(0.5^2))
p2slope2=(log10(errors(1,2))-log10(errors(1,3)))/(log10(0.5^2)-log10(0.5^3))
p2slope3=(log10(errors(1,3))-log10(errors(1,4)))/(log10(0.5^3)-log10(0.5^4))

p3slope1=(log10(errors(2,1))-log10(errors(2,2)))/(log10(0.5^1)-log10(0.5^2))
p3slope2=(log10(errors(2,2))-log10(errors(2,3)))/(log10(0.5^2)-log10(0.5^3))
p3slope3=(log10(errors(2,3))-log10(errors(2,4)))/(log10(0.5^3)-log10(0.5^4))

p4slope1=(log10(errors(3,1))-log10(errors(3,2)))/(log10(0.5^1)-log10(0.5^2))
p4slope2=(log10(errors(3,2))-log10(errors(3,3)))/(log10(0.5^2)-log10(0.5^3))
p4slope3=(log10(errors(3,3))-log10(errors(3,4)))/(log10(0.5^3)-log10(0.5^4))

% % p5slope1=(log10(errors(4,1))-log10(errors(4,2)))/(log10(0.5^1)-log10(0.5^2))
% % p5slope2=(log10(errors(4,2))-log10(errors(4,3)))/(log10(0.5^2)-log10(0.5^3))
% % p5slope3=(log10(errors(4,3))-log10(errors(4,4)))/(log10(0.5^3)-log10(0.5^4))

p2postslope1=(log10(errorspost(1,1))-log10(errorspost(1,2)))/(log10(0.5^1)-log10(0.5^2))
p2postslope2=(log10(errorspost(1,2))-log10(errorspost(1,3)))/(log10(0.5^2)-log10(0.5^3))
p2postslope3=(log10(errorspost(1,3))-log10(errorspost(1,4)))/(log10(0.5^3)-log10(0.5^4))

p3postslope1=(log10(errorspost(2,1))-log10(errorspost(2,2)))/(log10(0.5^1)-log10(0.5^2))
p3postslope2=(log10(errorspost(2,2))-log10(errorspost(2,3)))/(log10(0.5^2)-log10(0.5^3))
p3postslope3=(log10(errorspost(2,3))-log10(errorspost(2,4)))/(log10(0.5^3)-log10(0.5^4))

p4postslope1=(log10(errorspost(3,1))-log10(errorspost(3,2)))/(log10(0.5^1)-log10(0.5^2))
p4postslope2=(log10(errorspost(3,2))-log10(errorspost(3,3)))/(log10(0.5^2)-log10(0.5^3))
p4postslope3=(log10(errorspost(3,3))-log10(errorspost(3,4)))/(log10(0.5^3)-log10(0.5^4))

% % p5postslope1=(log10(errorspost(4,1))-log10(errorspost(4,2)))/(log10(0.5^1)-log10(0.5^2))
% % p5postslope2=(log10(errorspost(4,2))-log10(errorspost(4,3)))/(log10(0.5^2)-log10(0.5^3))
% % p5postslope3=(log10(errorspost(4,3))-log10(errorspost(4,4)))/(log10(0.5^3)-log10(0.5^4))







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



















