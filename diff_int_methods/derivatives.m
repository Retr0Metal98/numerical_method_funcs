function derivatives
%DERIVATIVES Has subfunctions for approximating 1st & 2nd derivatives
end

function df = get_df_by_fitting_func(x,y)
% input corresponds to n points
% if n = 3 -> fit quadratic y -> get linear dy/dx
% if n = 4 -> fit cubic y -> get quadratic dy/dx

n = length(x);
if n == 3
    % y = a*x^2 + b*x + c
    matrix = [x.^2 x x.*0+1];
    vector = y;
    R = rref_manual([matrix vector]);
    sol = R(:,end);
    a = sol(1); b = sol(2);
    % dy/dx = 2*a*x + b
    df = @(x) 2*a*x + b;
end
if n == 4
    % y = a*x^3 + b*x^2 + c*x + d
    matrix = [x.^3 x.^2 x x.*0+1];
    vector = y;
    R = rref_manual([matrix vector]);
    sol = R(:,end);
    a = sol(1); b = sol(2); c = sol(3);
    % dy/dx = 3*a*x^2 + 2*b*x + c
    df = @(x) 3*a*x^2 + 2*b*x + c;
end
end

function dft1 = first_2pt_backward(h,C0,C1)
% 2 point backward approximation for 1st derivative (used at end boundary)
% C0 = f(t0), C1 = f(t1), t0 = t1-h
dft1 = (C1-C0)/h;
end

function dft1 = first_2pt_forward(h,C1,C2)
% 2 point forward approximation for 1st derivative (used at start boundary)
% C1 = f(t1), C2 = f(t2), t2 = t1+h
dft1 = (C2-C1)/h;
end

function dft1 = first_central(h,C0,C2)
% central approximation for 1st derivative 
% C0 = f(t0), C2 = f(t2), t2 = t1+h, t0 = t1-h
dft1 = (C2-C0)/(2*h);
end

function dft2 = first_3pt_backward(h,C0,C1,C2)
% 3 point backward approximation for 1st derivative (used at end boundary)
% C0 = f(t0), C1 = f(t1), C2 = f(t2), t1 = t2-h, t0 = t2-2*h
dft2 = (C0-4*C1+3*C2)/(2*h);
end

function dft1 = first_3pt_forward(h,C1,C2,C3)
% 3 point forward approximation for 1st derivative (used at start boundary)
% C1 = f(t1), C2 = f(t2), C3 = f(t3), t2 = t1+h, t3 = t1+2*h
dft1 = (-3*C1+4*C2-C3)/(2*h);
end

function d2ft2 = second_backward(h,C0,C1,C2)
% 3 point backward approximation for 2nd derivative (used at end boundary)
% C0 = f(t0), C1 = f(t1), C2 = f(t2), t1 = t2-h, t0 = t2-2*h
d2ft2 = (C0-2*C1+C2)/(h^2);
end

function d2ft1 = second_forward(h,C1,C2,C3)
% 3 point forward approximation for 2nd derivative (used at start boundary)
% C1 = f(t1), C2 = f(t2), C3 = f(t3), t2 = t1+h, t3 = t1+2*h
d2ft1 = (C1-2*C2+C3)/(h^2);
end

function d2ft1 = second_central(h,C0,C1,C2)
% 3 point forward approximation for 2nd derivative (used at start boundary)
% C0 = f(t0), C1 = f(t1), C2 = f(t2), t0 = t1-h, t2 = t1+h
d2ft1 = (C0-2*C1+C2)/(h^2);
end

