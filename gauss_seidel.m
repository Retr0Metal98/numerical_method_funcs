function [x1, flag, enorm, err_list, err_count] = gauss_seidel(A,b,x0,num_iter)
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
zero_threshold = 1e-6;
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
