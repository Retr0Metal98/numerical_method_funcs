function NRG_test
v1 = [1 1]'; % initial values of input variables (x1, y1)
f = @(v) fun(v);
df = @(v) diff_fun(v);
num_iter = 1000;
% get the zeros in v2
[v2, flag, f2, curr_iter] = newton_raphson_general(v1,f,df,num_iter)
manual_example % for comparison with manually defining the functions individually
end

function manual_example
f = @(x,y) x^2 + x*y - 10;
g = @(x,y) y + 3*x*y^2 - 57;
 
dfdx = @(x,y) 2*x + y; dfdy = @(x,y) x;
dgdx = @(x,y) 3*y^2; dgdy = @(x,y) 1 + 3*x*2*y;
 
num = 1000; x1 = 1; y1 = 1;
 
for i=1:num
    A(1,1) = dfdx(x1,y1); A(1,2) = dfdy(x1,y1);
    A(2,1) = dgdx(x1,y1); A(2,2) = dgdy(x1,y1);
    b(1,1) = -f(x1,y1); b(2,1) = -g(x1,y1);
    s = A\b;
    x2 = x1 + s(1); y2 = y1 + s(2);
    if abs(f(x2,y2)) < 1e-10 && abs(g(x2,y2)) < 1e-10, break, end
    x1 = x2; y1 = y2;
end
x2, y2, f(x2,y2), g(x2,y2) 
end

function f = fun(v)
% Sample function 'fun'
% the inputs to f,g,etc. are collected in row/column vector v
% example: v(1) = x, v(2) = y, etc.
x = v(1); y = v(2);
% ensure that f is a column vector
f(1,1) = x^2 + x*y - 10;      % = f(x,y)
f(2,1) = y + 3*x*y^2 - 57;    % = g(x,y)
end

function df = diff_fun(v)
% Sample function 'diff_fun' - outputs the derivatives to functions in
% 'fun' as a matrix ('m' functions with 'n' variables = m x n matrix)
% the inputs to dfdx,dfdy,dgdx,dgdy, etc. are collected in row/column vector v
% example: v(1) = x, v(2) = y, etc.
x = v(1); y = v(2);

df = zeros(2,2);
% arranged such that columns = input variables, rows = functions
df(1,1) = 2*x + y; df(1,2) = x;
df(2,1) = 3*y^2; df(2,2) = 1 + 3*x*2*y;
end