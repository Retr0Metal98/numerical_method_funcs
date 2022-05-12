function [x1, flag, enorm] = gauss_seidel(A,b,x0,num_iter)
% Gauss-Seidel iterative algorithm to solve Ax=b
% for square matrices A
[m,n] = size(A);
Dinv = zeros(m,n);
C = A;
for i=1:n
    Dinv(i,i) = 1/A(i,i);
    C(i,i) = 0;
end

flag = 0; % whether solution is found
zero_threshold = 1e-10;
x = x0;
[x_m, x_n] = size(x);
for i=1:num_iter
    oldx = x; % x from previous iteration saved for later comparison
    for j=1:x_m
        x(j) = Dinv(j,j)*(b(j)-C(j,:)*x);
    end
    e = b - A*x; % length of error vector
    enorm = (e'*e)^0.5; % or sqrt(sum(e.^2));
    
    % normalised length of error vector (divide by length of b):
    % enorm = ((e'*e)^0.5)/((b'*b)^0.5); % or sqrt(sum(e.^2))/sqrt(sum(b.^2));
    % length of projection vector:
    % dx = x - oldx;
    % enorm = (dx'*dx)^0.5; % or sqrt(sum(dx.^2));
    if enorm < zero_threshold
        flag = 1;
        break;
    end
end
x1 = x;
end


function [x1, flag, enorm, err_list, err_count] = gauss_seidel_return_errs(A,b,x0,num_iter)
% Gauss-Seidel iterative algorithm to solve Ax=b
% for square matrices A
[m,n] = size(A);
Dinv = zeros(m,n);
C = A;
for i=1:n
    Dinv(i,i) = 1/A(i,i);
    C(i,i) = 0;
end

flag = 0; % whether solution is found
zero_threshold = 1e-10;
x = x0;
[x_m, x_n] = size(x);
err_list = zeros(num_iter,1);
err_count = 0;
for i=1:num_iter
    for j=1:x_m
        x(j) = Dinv(j,j)*(b(j)-C(j,:)*x);
    end
    e = b - A*x; % error vector
    enorm = (e'*e)^0.5;
    err_count = err_count + 1;
    err_list(err_count,1) = enorm;
    if enorm < zero_threshold
        flag = 1;
        break;
    end
end
x1 = x;
err_list = err_list(1:err_count,1);
end
