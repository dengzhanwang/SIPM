function [blk, data, At ] = multitask_read2(fname1,fname2,fname3,fname4,fname5,n0,p)
if nargin < 5
    p = 1;
end
data.trainsetall = [];
data.trainlabelall = [];
for i = 1:p
datatable = readtable(fname1);
data_matrix1 = table2array(datatable);
feature = data_matrix1(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix1(:,12)/norm(data_matrix1(:,12));



data.train_set1 = feature(1:n0,:);
data.train_label1 = label(1:n0);

data.test_set1 = feature(n0+1:n0+200,:);
data.test_label1 = label(n0+1:n0+200);
%%

datatable2 = readtable(fname2, 'Delimiter', ';');
data_matrix2 = table2array(datatable2);
feature = data_matrix2(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix2(:,12)/norm(data_matrix2(:,12));



data.train_set2 = feature(1:n0,:);
data.train_label2 = label(1:n0);

data.test_set2 = feature(n0+1:n0+200,:);
data.test_label2 = label(n0+1:n0+200);

%%

datatable3 = readtable(fname3);
data_matrix3 = table2array(datatable3);
feature = data_matrix3(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix3(:,9)/norm(data_matrix3(:,9));

data.train_set3 = feature(1:n0,:);
data.train_label3 = label(1:n0);

indextest = randperm(n0, 200);
data.test_set3 = feature(indextest,:);
data.test_label3 = label(indextest);

%%

datatable4 = readtable(fname4);
datatable4 = datatable4(:,3:15);
data_matrix4 = table2array(datatable4);

feature = data_matrix4(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix4(:,13)/norm(data_matrix4(:,13));

data.train_set4 = feature(1:n0,:);
data.train_label4 = label(1:n0);

data.test_set4 = feature(n0+1:n0+200,:);
data.test_label4 = label(2*n0+1:2*n0+200);

%%

% datatable5 = readtable(fname5);
datatable5 = importdata(fname5);
datatable5 =datatable5.data;
data_matrix5 = [ones(size( datatable5,1),1),datatable5 ];


feature = data_matrix5(:,1:8);
feature = feature./vecnorm(feature, 2, 1);
label = data_matrix2(:,9)/norm(data_matrix2(:,9));

data.train_set5 = feature(1:n0,:);
data.train_label5 = label(1:n0);

indextest = randperm(n0, 200);
data.test_set5 = feature(indextest,:);
data.test_label5 = label(indextest);

data.trainsetall = [data.trainsetall;data.train_set1; data.train_set2;data.train_set3;data.train_set4;data.train_set5 ];
data.trainlabelall = [data.trainlabelall;data.train_label1; data.train_label2;data.train_label3;data.train_label4;data.train_label5 ];
end


m = 5*p;

blk{1,1} = 's';
blk{1,2} = m;

Acell = cell(1, 1);

% tr(X) = K
Acell{1} = speye(m);

At = svec_sdpnal(blk, Acell, 1);

b = [1];

data.b = b;

end