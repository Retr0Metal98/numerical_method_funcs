function Null_A = null_manual(A)
%NULL_MANUAL returns the null space of A
A = rref_manual(A); % ensuring that A is in RREF before starting
[A_m,A_n] = size(A);
zero_threshold = 1e-10;

row = 1;
Fci = 1; % index of starting column of F

% keep track of column swaps, 
% these need to be applied to the final Null_A after adding the -I matrix
swap_count = 0;
swap_orig_cols = zeros(1,A_n);
swap_new_cols = zeros(1,A_n);
% locate the columns not part of identity
    for col = 1:A_n
        % if current col is not an identity col 
        % where the row'th element is a one, attempt swapping with
        % subsequent columns 
        if is_identity_col_for_row(A(:,col),row) == 0
            [A,orig_col,new_col] = swap_to_I_col(A,row,col);
            if orig_col ~= new_col
                swap_count = swap_count + 1;
                swap_orig_cols(1,swap_count) = orig_col;
                swap_new_cols(1,swap_count) = new_col;
            end
        end
        % if even after swapping, the current col is not the required
        % identity col, the identity part of the matrix is fully found
        % and subsequent cols form the 'F' part of the matrix
        if is_identity_col_for_row(A(:,col),row) == 0
            break;
        end
        Fci = Fci + 1;
        row = row + 1;
        if row > A_m
            break;
        end
    end 
F_cols = A(:,Fci:A_n);
[F_m, F_n] = size(F_cols);

% locate the end of 'F' part of null columns
non_zero_row = 0; % whether current row is has non-zeros
% check each row from bottom to top
    for nFr = F_m:-1:1
        for c = 1:F_n
            % if the current row has any non-zero values
            % then the last row of F has been reached
            if abs(F_cols(nFr,c)) > zero_threshold
                non_zero_row = 1;
            end
        end 
        if non_zero_row == 1
            break
        end 
    end
F = F_cols(1:nFr,:);

% create negative identity matrix of size (F_n,F_n)
minus_I = [-1];
    for k = 2:F_n
        minus_I(k,k) = -1;
    end 
% concatenate F with the -ve identity matrix (-I) vertically
Null_A = [F;minus_I];
    if swap_count > 0
        % if any columns were swapped to obtain F
        % the corresponding rows of Null_A must be swapped
        % to obtain the correct null space vectors
        % e.g.: if columns 2 & 3 were swapped, rows 2 & 3 must be swapped
        for q = 1:swap_count
            orig_idx = swap_orig_cols(1,q);
            new_idx = swap_new_cols(1,q);
            Null_A = swap_any_rows(Null_A,orig_idx,new_idx);
        end
    end
end

function [A,orig_col,swap_col] = swap_to_I_col(A,piv_row,col)
% swap the given col of A with a col to the right which is an
% identity column whose piv_row is 1.
[m,n] = size(A);
left = A(:,1:col-1); % keep columns to the left of target column unchanged
right = A(:,col:n); % columns to the right of target column
[right_m, right_n] = size(right);
% keep track of which columns of A are swapped
% to start, no swap has occured
orig_col = col;
swap_col = col;
    for c = 1:right_n
        % go through each column to the right of target column
        % if the column is an identity column with the correct row = 1
        % then perform the swap
        if is_identity_col_for_row(right(:,c),piv_row) == 1
            right = swap_any_cols(right,1,c);
            % column which is swapped with original column of A
            % has index c with respect to right
            % but need to add n-right_n to get index with respect to A
            swap_col = c+(n-right_n);
            break;
        end
    end 
A = [left right];
end

function res = is_identity_col_for_row(col_vec,target_row)
% function to check if the given col an identity column
% where the pivot (=1) is in the target_row
res = 1;
zero_threshold = 1e-10;
[num_rows,num_cols] = size(col_vec);
    for r = 1:num_rows
        % check each row of the column vector
        if r == target_row && col_vec(r,1) ~= 1
            % if the row index matches the target row (example: 2nd row)
            % and the value ~= 1, then the column is not the required
            % identity column
            res = 0;
            break;
        end 
        if r ~= target_row && abs(col_vec(r,1)) > zero_threshold
            % if any of the other rows is not zero, then the column is not
            % the required identity column
            res = 0;
            break;
        end 
    end
end

function A = swap_any_rows(A,row_a_idx,row_b_idx)
% function to swap any two rows of a matrix
temp = A(row_a_idx,:);
A(row_a_idx,:) = A(row_b_idx,:);
A(row_b_idx,:) = temp;
end

function A = swap_any_cols(A,col_a_idx,col_b_idx)
% function to swap any two columns of a matrix
temp = A(:,col_a_idx);
A(:,col_a_idx) = A(:,col_b_idx);
A(:,col_b_idx) = temp;
end
