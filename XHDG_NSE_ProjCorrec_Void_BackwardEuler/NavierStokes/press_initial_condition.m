function p_i = press_initial_condition(X,time)


x = X(:,1);
y = X(:,2);
C = 0;

%% here constant could be time or a constant.

%% initial condition for pressure:
%p_i = x+y; %p(x,y,time) = x+y+time, for time=0, p(x,y,time=0)=x+y.
%
%% Benchmark example from Ueckermann & Lermusiaux
p_i = C.*(x+y);
%p_i = sin(C).*cos(pi.*x).*sin(pi.*y);


%p_i = -x;

%disp('Hola')



