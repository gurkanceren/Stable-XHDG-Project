function printSlopes(x,y,eps)

n = length(x);
for i=1:n-1
    a=x(i); b=x(i+1); fa=y(i); fb=y(i+1);
    L = b-a; H = fb-fa; slope = abs(H/L);
    xt = a+L/2; yt = fa+H*(0.5+eps);
    text(xt,yt,sprintf('%0.2g',slope));
end