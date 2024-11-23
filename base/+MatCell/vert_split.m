%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 13:54:20
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 09:48:33
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = vert_split(A, block_size)
    % vert_split - Split rows of the matrix A into blocks
    % Input:
    %   A - matrix to be split
    %   block_size - array containing the number of rows for each block
    %
    % Output:
    %   out - cell containing the split blocks

    % Check if the sum of block_size equals the number of rows in A
    assert(sum(block_size) == size(A, 1), 'The sum of block_size must equal the number of rows in A');

    % Initialize the MatCell outect with the length of block_size
    out = cell(length(block_size), 1);
    if length(block_size) == 1
        out.data = A;
        return
    end
    % Split the rows of A into blocks and store them in the MatCell outect
    row_start = 1;
    for k = 1:length(block_size)
        row_end = row_start + block_size(k) - 1;
        out{k} = A(row_start:row_end, :);
        row_start = row_end + 1;
    end
    
end