function [x3,flag,y3,curr_iter] = secant(x1,x2,f,num_iter)
%SECANT
zero_threshold = 1e-10;
flag = 0;

for curr_iter=1:num_iter
    matrix = [x1 1
              x2 1];
    vector = [f(x1)
              f(x2)];
    s = matrix\vector;
    a = s(1); b = s(2);
    x3 = -b/a;
    y3 = f(x3);
    if abs(y3) < zero_threshold
        flag = 1;
        break;
    end 
    x1 = x2; x2 = x3;
end
end

