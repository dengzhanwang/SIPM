
function B = rearrangeColumns(A)
    % 确保输入矩阵是 5x7
    if size(A,1) ~= 5 || size(A,2) ~= 7
        error('Input matrix must be a 5x7 matrix.');
    end
    
    % 创建一个新矩阵，初始与 A 相同
    B = A;
    
    % 把第五列移到第二列
    B(:,2) = A(:,5);
    
    % 把原来的第二、三、四列依次向后移动
    B(:,3) = A(:,2);
    B(:,4) = A(:,3);
    B(:,5) = A(:,4);
    
    % 其他列保持不变
    B(:,1) = A(:,1);
    B(:,6) = A(:,6);
    B(:,7) = A(:,7);
end
