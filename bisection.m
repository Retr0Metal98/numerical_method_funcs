function [x3,flag,y3,curr_iter] = bisection(x1,x2,f)
%BISECTION
zero_threshold = 1e-10;
num_iter = 1e6;
flag = 0;
% run the iterations of the bisection algorithm
for curr_iter=1:num_iter
    x3 = (x1+x2)/2;
    y3 = f(x3);
    % stop iteration if |f(x)| close enough to zero
    if abs(y3) < zero_threshold
        flag = 1;
        break;
    end
    % proceed to next iteration otherwise
    % new x1 & x2 such that f(x1) & f(x2) have opposite signs
    if f(x3)*f(x1) < 0
        x2 = x3;
    else
        x1 = x3;
    end
end 
end

