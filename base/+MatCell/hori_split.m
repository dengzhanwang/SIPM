%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 13:52:38
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 09:47:34
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = hori_split(A, block_size)
    % hori_split - Split columns of the matrix A into blocks
    % Input:
    %   A - matrix to be split
    %   block_size - array containing the number of columns for each block
    %
    % Output:
    %   out - cell containing the split blocks

    % Check if the sum of block_size equals the number of columns in A
    assert(sum(block_size) == size(A, 2), 'The sum of block_size must equal the number of columns in A');

    % Initialize the MatCell outect with the length of block_size
    out = cell(length(block_size), 1);
    if length(block_size) == 1
        out.data = A;
        return
    end
    % Split the columns of A into blocks and store them in the MatCell outect
    col_start = 1;
    for k = 1:length(block_size)
        col_end = col_start + block_size(k) - 1;
        out{k} = A(:, col_start:col_end);
        col_start = col_end + 1;
    end
end