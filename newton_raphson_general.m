function [x2,flag,f2,curr_iter] = newton_raphson_general(x1,f,df,num_iter)
%NEWTON_RAPHSON for any number of functions & variables: see NRG_test.m for an example on how to use this function
zero_threshold = 1e-10;
flag = 0;
for curr_iter=1:num_iter
    A = df(x1);
    b = -f(x1);
    % solve Ax = b, where x = col vector: each entry = x2 - x1
    Ab_aug_R = rref_manual([A b]);
    s = Ab_aug_R(:,end); % get last column
    x2 = x1+s; % update 
    f2 = f(x2);
    if abs(f2)<zero_threshold
        flag = 1;
        break;
    end
    x1 = x2;
end
end