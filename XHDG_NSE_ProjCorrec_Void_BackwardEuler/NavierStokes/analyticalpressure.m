function p = analyticalpressure(X,time)


x = X(:,1);
y = X(:,2);
C = 0;

%% here constant could be time or a constant.

%% initial condition for pressure:
%p_i = x+y; %p(x,y,time) = x+y+time, for time=0, p(x,y,time=0)=x+y.
%p = 0+x*0; %x+y+C;
p = 1.*(x+y+C);
%
%% Benchmark example from Ueckermann & Lermusiaux
%p = sin(time).*cos(pi.*x).*sin(pi.*y);
%
%
%p_i = -x;

%disp('Hola')



