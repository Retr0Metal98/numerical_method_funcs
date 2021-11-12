function total_area = trapezoid(x,y)
%TRAPEZOID 
% x,y are two column vectors assumed to be same size
% check with built-in function `trapz(t,y)`
total_area = 0;
[m,n] = size(x);
for j=1:m-1
    e1 = y(j);
    e2 = y(j+1);
    h = x(j+1)-x(j);
    sub_area = ((e1+e2)/2)*h;
    total_area = total_area + sub_area;
end
end

