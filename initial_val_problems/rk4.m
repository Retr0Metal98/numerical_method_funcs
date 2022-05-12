function [func_data] = rk4(df_fun,t,init_vals)
%RK4
% Higher order ODE -> convert to system of 1st order ODEs
k = length(t);
h = t(2) - t(1);
v = size(init_vals,2);
func_data = zeros(k,v); 
% k rows for k indep var values, v cols for v dependent vars
func_data(1,:) = init_vals;

for i=1:k-1
    t1 = t(i); f1 = func_data(i,:);
    % estimate of 2nd point using 1st point's derivative to get area
    diff1 = df_fun(t1,f1);
    f2_1 = f1 + diff1.*h./2;
    t2 = t1 + h/2;
    diff2_1 = df_fun(t2,f2_1);
    % re-estimate 2nd point using previous estimate of 2nd point to get area
    f2_2 = f1 + diff2_1.*h./2;
    diff2_2 = df_fun(t2,f2_2);
    % use 2nd estimate of 2nd point to get 3rd point
    f3 = f1 + diff2_2.*h;
    t3 = t1 + h;
    diff3 = df_fun(t3,f3);
    diff2 = (diff2_1 + diff2_2)./2;
    % use 3 derivative values to get true f using 1/3 Simpson's rule
    f3 = f1 + h./6.*(diff1+4.*diff2+diff3);
    func_data(i+1,:) = f3;
end
end

% Sample function 'fun'
function dfdt = fun(t,f)
dfdt(1) = 0.66*(1-f(1))*(1-f(2));
dfdt(2) = 2*(1-f(1))*(1-f(2));
dfdt(3) = -(1-f(1))*f(3);
end