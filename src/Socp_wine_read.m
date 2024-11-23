function [blk, data, At ] = Socp_wine_read(fname)
% data_red = readtable('fname', 'Delimiter', ';');
% data_matrix = table2array(data_red);


datatable = readtable(fname, 'Delimiter', ';');
data_matrix = table2array(datatable);
feature = data_matrix(:,1:11);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix(:,12);
Sigma = (cov(feature,1))^(1/2)*10000;

data.train_set = feature(1:500,:);
data.train_label = label(1:500);
data.test_set = feature(501:end,:);
data.test_label = label(501:end);

data.train_miss = generate_random_0_1_matrix(250,size(data.train_set,2));


% data.train_miss = rand(250,size(data.train_set,2)) > 0.2;
data.train_miss = [data.train_set(1:250,:); data.train_miss.*data.train_set(251:end,:)];
data.train_miss = fill_zeros_with_mean(data.train_miss); 


n = size(feature,2);
At{1} = [zeros(n,1),eye(n),zeros(n,1),-Sigma]';
blk{1,1} = 'q';
blk{1,2} = n;

blk{2,1} = 'q';
blk{2,2} = n;



end
function matrix = generate_random_0_1_matrix(rows, cols)
    matrix = ones(rows, cols);  % 创建全1矩阵
    for i = 1:rows
        zero_indices = randperm(cols, 2);  % 随机选择两个位置作为0的位置
        matrix(i, zero_indices) = 0;
    end
end

function result = fill_zeros_with_mean(matrix)
    % 找到非零元素
    non_zero_elements = matrix(matrix ~= 0);
    
    % 计算非零元素的均值
    mean_value = mean(non_zero_elements);
    
    % 将矩阵中的0替换为均值
    result = matrix;
    result(result == 0) = mean_value;
end