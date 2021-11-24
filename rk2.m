function [func_data] = rk2(df_fun,t,init_vals)
%RK2
% Higher order ODE -> convert to system of 1st order ODEs
k = length(t);
h = t(2) - t(1);
v = size(init_vals,2);
func_data = zeros(k,v); 
% k rows for k indep var values, v cols for v dependent vars
func_data(1,:) = init_vals;

for i=1:k-1
    t1 = t(i); f1 = func_data(i,:);
    diff1 = df_fun(t1,f1);
    f2 = f1 + diff1.*h; t2 = t1 + h;
    diff2 = df_fun(t2,f2);
    f2 = f1 + (diff1 + diff2)./2*h;
    func_data(i+1,:) = f2;
end
end

% Sample function 'fun'
function dfdt = fun(t,f)
dfdt(1) = f(1)*t^2;
%dfdt(2) = 2*(1-f(1))*(1-f(2));
%dfdt(3) = -(1-f(1))*f(3);
end