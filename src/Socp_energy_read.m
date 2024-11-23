function [blk, data, At ] = Socp_energy_read(fname,n0)
try
datatable = readtable(fname);
data_matrix = table2array(datatable);
catch
data_matrix = importdata(fname);
end
feature = data_matrix(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix(:,end);
Sigma = (cov(feature,1)+1e-8*eye(8))^(1/2)*10000  ;


data.train_set = feature(1:n0,:);
data.train_label = label(1:n0);
data.test_set = feature(n0+1:end,:);
data.test_label = label(n0+1:end);
data.train_miss = generate_random_0_1_matrix(ceil(n0/2),size(data.train_set,2));
data.train_miss = [data.train_set(1:ceil(n0/2),:); data.train_miss.*data.train_set(ceil(n0/2)+1:end,:)];
data.train_miss = fill_zeros_with_mean(data.train_miss); 




n = size(feature,2);
At{1} = [zeros(n,1),eye(n),zeros(n,1),-Sigma]';
blk{1,1} = 'q';
blk{1,2} = n;

blk{2,1} = 'q';
blk{2,2} = n;



end

function matrix = generate_random_0_1_matrix(rows, cols)
    matrix = ones(rows, cols); 
    for i = 1:rows
        zero_indices = randperm(cols, 2); 
        matrix(i, zero_indices) = 0;
    end
end

function result = fill_zeros_with_mean(matrix)

    non_zero_elements = matrix(matrix ~= 0);
    

    mean_value = mean(non_zero_elements);
    

    result = matrix;
    result(result == 0) = mean_value;
end