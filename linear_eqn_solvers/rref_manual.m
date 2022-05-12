function A = rref_manual(A)
%RREF_MANUAL returns RREF of matrix input A
piv_row = 1;
[m,n] = size(A);
zero_threshold = 1e-10;
    for piv_col = 1:n
        pivot = A(piv_row,piv_col);
        if abs(pivot) < zero_threshold
            A = swap_to_piv_row(A,piv_row,piv_col); % swap row w/ zero pivot with other row
            pivot = A(piv_row,piv_col); % update pivot
        end
        if abs(pivot) < zero_threshold
            % if pivot is still zero, move to next column, 
            % but stay on same row
            continue 
        end
        % normalize the pivot row by pivot value 
        % divide throughout to make pivot value = 1
        A(piv_row,:) = A(piv_row,:)/pivot;
        for r=1:m
            % for all rows other than pivot row
            % make the leading entry zero by subtracting the product of
            % leading entry and normalized pivot row 
            if r ~= piv_row
                A(r,:) = A(r,:) - A(r,piv_col)*A(piv_row,:);
            end 
        end
        if piv_row == m
            % if the last row has already been processed
            % RREF algo is completed
            break
        end
        % move to next row after fixing the pivot for this column
        % i.e., move diagonally down to next row
        piv_row = piv_row+1;
    end
end

function A = swap_to_piv_row(A,piv_row,piv_col)
[m,n] = size(A);
upper = A(1:piv_row-1,:); % keep rows above pivot unchanged
lower = A(piv_row:m,:); % pivot row + rows below it
piv_col_abs_vals = abs(lower(:,piv_col));
[max_abs_val,max_val_idx] = max(piv_col_abs_vals);
% swap the row of the pivot element (1st row of lower) with row below pivot 
% which has max absolute value among other rows below pivot
lower = swap_any_rows(lower,1,max_val_idx);
% recombine the upper rows (unchanged) with the changed lower rows
A = [upper;lower];
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