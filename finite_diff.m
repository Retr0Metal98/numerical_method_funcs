function [y,yref] = finite_diff
%FINITE_DIFF Example of finite difference method
% System of ODEs -> combine into one higher order ODE (only linear ODEs considered)

% ODE: y'' + p(x)*y' + q(x)*y = r(x)
% example: y'' - 2*y'+ 3*y = 5, y'(0) = 1, y'(1) = 0.
h = 0.1;
x = 0:h:1;
k = length(x);
A = zeros(k,k); % coeff matrix
b = zeros(k,1); % RHS of eqns

% Boundary conditions - see boundary_conditions subfunction for examples.
% here y'(0) = 1, y'(1) = 0.
A(1,1) = -3/(2*h); A(1,2) = 4/(2*h); A(1,3) = -1/(2*h); b(1,1) = 1;
A(k,k-2) = 1/(2*h); A(k,k-1) = -4/(2*h); A(k,k) = 3/(2*h); b(k,1) = 0;

for i=2:k-1
    p = -2; q = 3; r = 5; % p, q, r are functions of x
    % using central method for y' & y''
    % y'' = (y(i+1)-2*y(i)+y(i+1)/h^2
    % p*y' = p*(y(i+1)-y(i-1))/(2*h)
    % q*y  = q*y(i)
    % rearranging to isolate y(i-1), y(i), and y(i+1):
    A(i,i-1) = 1/h^2 - p/(2*h); % = y(i-1)
    A(i,i) = -2/h^2 + q;        % = y(i)
    A(i,i+1) = 1/h^2 + p/(2*h); % = y(i+1)
    b(i,1) = r;
end

% get y by solving Ay = b, using rref_manual, jacobi, gauss-seidel, etc.
Ab_aug_R = rref_manual([A b]);
y = Ab_aug_R(:,end);
yref = A\b;

end

function boundary_conditions
A = zeros(k,k); % coeff matrix
b = zeros(k,1); % RHS of eqns

u = 3; v = 7;
% for 1st BC: use 1st row, for 2nd BC: use kth row
% for y = u (BC1), y = v (BC2):
A(1,1) = 1; b(1,1) = u;
A(k,k) = 1; b(k,1) = v;

% for y' = u (BC1) & y' = v (BC2)
% use 3 point forward for BC1: y' = (-3*y(1)+4*y(2)-y(3))/(2*h)
% use 3 point backward for BC2: y' = (y(k-2)-4*y(k-1)+3*y(k))/(2*h)
A(1,1) = -3/(2*h); A(1,2) = 4/(2*h); A(1,3) = -1/(2*h); b(1,1) = u;
A(k,k-2) = 1/(2*h); A(k,k-1) = -4/(2*h); A(k,k) = 3/(2*h); b(k,1) = v;

% for y'' = u (BC1) & y'' = v (BC2)
% use 2nd deriv. forward for BC1: y'' = (y(1)-2*y(2)+y(3))/(h^2)
% use 2nd deriv. backward for BC2: y'' = (y(k-2)-2*y(k-1)+y(k))/(h^2)
A(1,1) = 1/h^2; A(1,2) = -2/h^2; A(1,3) = 1/h^2; b(1,1) = u;
A(k,k-2) = 1/h^2; A(k,k-1) = -2/h^2; A(k,k) = 1/h^2; b(k,1) = v;
% can mix & match the above as needed
end
