function [blk, data, At ] = multitask_read(fname1,fname2,fname3)

datatable = readtable(fname1);
data_matrix1 = table2array(datatable);
feature = data_matrix1(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix1(:,12);

data.train_set1 = feature(1:500,:);
data.train_label1 = label(1:500);

datatable2 = readtable(fname2, 'Delimiter', ';');
data_matrix2 = table2array(datatable2);
feature = data_matrix2(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix2(:,12);

data.train_set2 = feature(1:500,:);
data.train_label2 = label(1:500);

datatable3 = readtable(fname3);
data_matrix3 = table2array(datatable3);
feature = data_matrix3(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix3(:,9);

data.train_set3 = feature(1:500,:);
data.train_label3 = label(1:500);



feature = data_matrix1(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix1(:,12);

data.train_set4 = feature(501:1000,:);
data.train_label4 = label(501:1000);


feature = data_matrix2(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix1(:,12);

data.train_set5 = feature(501:1000,:);
data.train_label5 = label(501:1000);


m = 5;

blk{1,1} = 's';
blk{1,2} = m;

Acell = cell(1, 1);

% tr(X) = K
Acell{1} = speye(m);

At = svec_sdpnal(blk, Acell, 1);

b = [1];

data.b = b;

end