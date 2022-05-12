function [x2,flag,y2,curr_iter] = newton_raphson(x1,f,df,num_iter)
%NEWTON_RAPHSON 
zero_threshold = 1e-10;
flag = 0;
for curr_iter=1:num_iter
    if isa(df,'function_handle') ~= 1
        % if df is not a function handle, use central difference
        % approximation to find derivative of function at x1
        h = 1e-4;
        df_x1 = (f(x1+h)-f(x1-h))/(2*h);
    else
        df_x1 = df(x1);
    end
    x2 = x1 - f(x1)/df_x1;
    y2 = f(x2);
    if abs(y2)<zero_threshold
        flag = 1;
        break;
    end
    x1 = x2;
end
end

