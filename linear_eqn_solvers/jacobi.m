function [x1, flag, enorm] = jacobi(A,b,x0,num_iter)
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
for i=1:num_iter
    x1 = Dinv*(b-C*x0); % new x
    e = b - A*x1; % length of error vector
    enorm = (e'*e)^0.5; % or sqrt(sum(e.^2));
    
    % normalised length of error vector (divide by length of b):
    % enorm = ((e'*e)^0.5)/((b'*b)^0.5); % or sqrt(sum(e.^2))/sqrt(sum(b.^2));
    % length of projection vector:
    % dx = x1 - x0;
    % enorm = (dx'*dx)^0.5; or sqrt(sum(dx.^2));
    if enorm < zero_threshold 
        flag = 1;
        break;
    end
    x0 = x1; % update old x to new x
end
end

function [x1, flag, enorm, err_list, err_count] = jacobi_return_errs(A,b,x0,num_iter)
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
    if enorm < zero_threshold % length of error vector
        flag = 1;
        break;
    end
    x0 = x1; % update old x to new x
end
err_list = err_list(1:err_count,1);
end
