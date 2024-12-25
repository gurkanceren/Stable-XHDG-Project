load Laplace/errors_P3_tau1
errorHDG_tau1=errorHDG; errorHDGpost_tau1=errorHDGpost;

load Laplace/errors_P3_tau100
errorHDG_tau100=errorHDG; errorHDGpost_tau100=errorHDGpost;

load Laplace/errors_P3_tau1singleface
errorHDG_tau1sf=errorHDG; errorHDGpost_tau1sf=errorHDGpost;

load Laplace/errors_P3_tau100singleface
errorHDG_tau100sf=errorHDG; errorHDGpost_tau100sf=errorHDGpost;

h=0.5*0.5.^[0:4];

figure(11),clf
plot(log10(h),log10(errorHDG_tau1),'bo-',log10(h),log10(errorHDGpost_tau1),'b-*'...
    ,log10(h),log10(errorHDG_tau100),'ms-',log10(h),log10(errorHDGpost_tau100),'m-x')
xlabel('log_{10}(h)'), ylabel('log_{10}(L2 error)')
legend('HDG \tau=1','HDG post \tau=1','HDG \tau=100','HDG post \tau=100',4)
title('HDG with \tau in all faces')

figure(12),clf
plot(log10(h),log10(errorHDG_tau1sf),'bo-',log10(h),log10(errorHDGpost_tau1sf),'b-*'...
    ,log10(h),log10(errorHDG_tau100sf),'ms-',log10(h),log10(errorHDGpost_tau100sf),'m-x')
xlabel('log_{10}(h)'), ylabel('log_{10}(L2 error)')
legend('HDG \tau=1','HDG post \tau=1','HDG \tau=100','HDG post \tau=100',4)
title('HDG with \tau in one face per element (SINGLE FACE)')

figure(13),clf
loglog(h,errorHDG_tau1sf,'bo-',h,errorHDGpost_tau1sf,'m-*','LineWidth',2)
xlabel('Element size'), ylabel('L2 error')
% legend('HDG \tau=1','HDG post \tau=1','HDG \tau=100','HDG post \tau=100',4)
% title('HDG with \tau in one face per element (SINGLE FACE)')

