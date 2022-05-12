function [total_area] = simpsons(t,y)
%SIMPSONS Numerical integration by Simpson's 1/3 & 3/8 rules
%   't' needs to be evenly distributed: use start:step:stop
k = length(t); % number of steps
h = t(2) - t(1); % step size
total_area = 0;
for i=1:k
    k2 = length(t);
    if k2 == 4 % use 3/8 Simpson's rule if 4 points left (occurs when number of steps is even)
        h = (t(4) - t(1))/2;
        area_slice = h/4 * (y(1)+3*y(2)+3*y(3)+y(4));
        total_area = total_area + area_slice;
        break
    end
    % otherwise use 1/3 Simpson's rule
    area_slice = h/3 * (y(1)+4*y(2)+y(3));
    total_area = total_area + area_slice;
    if k2 == 3, break, end
    % remove first 2 points and move to next iteration starting from 3rd point
    t = t(3:k2);
    y = y(3:k2);
end
end

