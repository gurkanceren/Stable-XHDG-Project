function u_i = initial_condition(X)


x = X(:,1);
y = X(:,2);


%% initial condition
u_i = x;
%u_i = x+y;
%u_i = x.^2+y.^2;
%u_i = x.*y;
%u_i = (x.^2).*y;
%u_i = x.*(y.^2);
%u_i = x.^3+y.^3;
%u_i = x.^4+y.^4;
%u_i = x.^5+y.^5;
%u_i = x.^6+y.^6;

%alpha_x = 0.01; alpha_y = 0.01;
%x_0 = 0.5; y_0 = 0.5;
%aa = (x-x_0).^2; bb = (y-y_0).^2;

%a = -aa/alpha_x;
%b = -bb/alpha_y;

%u_i = exp(a+b);


%disp('Hola')



