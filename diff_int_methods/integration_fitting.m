

function [I,s] = integration_fitting(x,y)
% input corresponds to n points
% if n = 2 -> fit linear y -> get quadratic integral function
% if n = 3 -> fit quadratic y -> get cubic integral function
% if n = 4 -> fit cubic y -> get 4th order integral function

n = length(x);
if n == 2
    % y = a*x + b % same as trapezoidal rule
    matrix = [x x.*0+1];
    vector = y;
    R = rref_manual([matrix vector]);
    sol = R(:,end);
    a = sol(1); b = sol(2); 
    % Y = a*x^2/2 + b*x
    I = @(x) a/2*x^2 + b*x;
    s = I(x(k)) - I(x(1));
end
if n == 3
    % y = a*x^2 + b*x + c
    matrix = [x.^2 x x.*0+1];
    vector = y;
    R = rref_manual([matrix vector]);
    sol = R(:,end);
    a = sol(1); b = sol(2); c = sol(3);
    % Y = a*x^3/3 + b*x^2/2 + c*x
    I = @(x) a/3*x^3 + b/2*x^2 + c*x;
    s = I(x(k)) - I(x(1));
end
if n == 4
    % y = a*x^3 + b*x^2 + c*x + d
    matrix = [x.^3 x.^2 x x.*0+1];
    vector = y;
    R = rref_manual([matrix vector]);
    sol = R(:,end);
    a = sol(1); b = sol(2); c = sol(3); d = sol(4);
    % Y = a*x^4/4 + b*x^3/3 + c*x^2/2 + d*x
    I = @(x) a/4*x^4 + b/3*x^3 + c/2*x^2 + d*x;
    s = I(x(k)) - I(x(1));
end
end