function [x1, flag, enorm, err_list, err_count] = jacobi(A,b,x0,num_iter)
% Jacobi iterative algorithm to solve Ax=b
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
err_list = zeros(num_iter,1);
err_count = 0;
for i=1:num_iter
    x1 = Dinv*(b-C*x0); % new x
    e = b - A*x1; % error vector
    enorm = (e'*e)^0.5;
    err_count = err_count + 1;
    err_list(err_count,1) = enorm;
    if enorm < zero_threshold
        flag = 1;
        break;
    end
    x0 = x1; % update old x to new x
end
err_list = err_list(1:err_count,1);
end
