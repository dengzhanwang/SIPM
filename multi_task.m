clear
close all;
test_dir = './';
seed = 2024;
rng(seed);
addpath(genpath(test_dir));
fvaltable = zeros(2,2,7);
errortable = zeros(2,2,7);
traintable = zeros(2,7,5);
testtable = zeros(2,7,5);
pall = [1 2];
n0all = [200 500];

p = 2;
for p = 1
    for k = 2
n0 = n0all(k);
[blk, data, At ] = multitask_read2('winequality-red.csv','winequality-white.csv','ENB2012_data.xlsx','AirQualityUCI.xlsx','abalone.data',n0,p);


opts.lambda1 = 1e-2;
opts.lambda2 = 1e-2;
opts.batchsize = n0;

opts.maxiter = 10000;
Amap = @(x) At'*x;
ATmap = @(y) At*y;
opts.Amap = Amap;
opts.ATmap = ATmap;
opts.n = 8;
opts.gamma = 0.1;
x0 = zeros(opts.n,p*5);
opts.Amap = @(X) AXmap(X, blk, At);
opts.ATmap = @(y) Atymap(y, blk, At);
b = data.b;

opts.numall =5*p;
opts.m = 5*p;
Omega0 = 1/opts.m*eye(opts.m,opts.m);
opts.Omega0 = Omega0;
opts.methods = 'sto';
opts.lr = @(i) min(max(3/i^(1/1.5),0.01),1);
opts.options = 1;


out1 = sto_ipm_multi_task3(blk,x0,data,At,b,opts);
out1.ferror = out1.fval(end)/abs(out1.fval(1));
fvaltable(p,k,1) = out1.ferror;
errortable(p,k,1) = out1.error(end)/max(out1.error);
traintable(k,1,1) = norm(data.train_set1*out1.W(:,1) - data.train_label1)/norm(data.train_label1);
testtable(k,1,1) = norm(data.test_set1*out1.W(:,1) - data.test_label1)/norm(data.test_label1);
traintable(k,1,2) = norm(data.train_set2*out1.W(:,2) - data.train_label2)/norm(data.train_label2);
testtable(k,1,2) = norm(data.test_set2*out1.W(:,2) - data.test_label2)/norm(data.test_label2);
traintable(k,1,3) = norm(data.train_set3*out1.W(:,3) - data.train_label3)/norm(data.train_label3);
testtable(k,1,3) =  norm(data.test_set3*out1.W(:,3) - data.test_label3)/norm(data.test_label3);
traintable(k,1,4) = norm(data.train_set4*out1.W(:,4) - data.train_label4)/norm(data.train_label4);
testtable(k,1,4) = norm(data.test_set4*out1.W(:,4) - data.test_label4)/norm(data.test_label4);
traintable(k,1,5) = norm(data.train_set5*out1.W(:,5) - data.train_label5)/norm(data.train_label5);
testtable(k,1,5) = norm(data.test_set5*out1.W(:,5) - data.test_label5)/norm(data.test_label5);

norm(out1.W'*out1.W/trace(out1.W'*out1.W) - out1.Omega{1},'fro')


opts.methods = 'mom';
opts.lr = @(i) min(max(4/i^(1/1.6),0.01),1);
out2 = sto_ipm_multi_task3(blk,x0,data,At,b,opts);
out2.ferror = out2.fval(end)/abs(out1.fval(1));

fvaltable(p,k,2) = out2.ferror;
errortable(p,k,2) = out2.error(end)/out2.error(1);

traintable(k,2,1) = norm(data.train_set1*out2.W(:,1) - data.train_label1)/norm(data.train_label1);
testtable(k,2,1) = norm(data.test_set1*out2.W(:,1) - data.test_label1)/norm(data.test_label1);
traintable(k,2,2) = norm(data.train_set2*out2.W(:,2) - data.train_label2)/norm(data.train_label2);
testtable(k,2,2) = norm(data.test_set2*out2.W(:,2) - data.test_label2)/norm(data.test_label2);
traintable(k,2,3) = norm(data.train_set3*out2.W(:,3) - data.train_label3)/norm(data.train_label3);
testtable(k,2,3) =  norm(data.test_set3*out2.W(:,3) - data.test_label3)/norm(data.test_label3);
traintable(k,2,4) = norm(data.train_set4*out2.W(:,4) - data.train_label4)/norm(data.train_label4);
testtable(k,2,4) = norm(data.test_set4*out2.W(:,4) - data.test_label4)/norm(data.test_label4);
traintable(k,2,5) = norm(data.train_set5*out2.W(:,5) - data.train_label5)/norm(data.train_label5);
testtable(k,2,5) = norm(data.test_set5*out2.W(:,5) - data.test_label5)/norm(data.test_label5);



opts.gamma = 0.9;
opts.methods = 'igt';
opts.lr = @(i) min(max(5/i^(1/1.6),0.01),1);

out3 = sto_ipm_multi_task3(blk,x0,data,At,b,opts);
out3.ferror = out3.fval(end)/abs(out1.fval(1));

fvaltable(p,k,3) = out3.ferror;
errortable(p,k,3) = out3.error(end)/out1.error(1);

traintable(k,3,1) = norm(data.train_set1*out3.W(:,1) - data.train_label1)/norm(data.train_label1);
testtable(k,3,1) = norm(data.test_set1*out3.W(:,1) - data.test_label1)/norm(data.test_label1);
traintable(k,3,2) = norm(data.train_set2*out3.W(:,2) - data.train_label2)/norm(data.train_label2);
testtable(k,3,2) = norm(data.test_set2*out3.W(:,2) - data.test_label2)/norm(data.test_label2);
traintable(k,3,3) = norm(data.train_set3*out3.W(:,3) - data.train_label3)/norm(data.train_label3);
testtable(k,3,3) =  norm(data.test_set3*out3.W(:,3) - data.test_label3)/norm(data.test_label3);
traintable(k,3,4) = norm(data.train_set4*out3.W(:,4) - data.train_label4)/norm(data.train_label4);
testtable(k,3,4) = norm(data.test_set4*out3.W(:,4) - data.test_label4)/norm(data.test_label4);
traintable(k,3,5) = norm(data.train_set5*out3.W(:,5) - data.train_label5)/norm(data.train_label5);
testtable(k,3,5) = norm(data.test_set5*out3.W(:,5) - data.test_label5)/norm(data.test_label5);


hold on

opts.radius = 400000;
opts.methods = 'recursiv';
opts.lr = @(i) min(max(14/i^(1/1.5),0.01),1);

out4 = sto_ipm_multi_task3(blk,x0,data,At,b,opts);
out4.ferror = out4.fval(end)/abs(out1.fval(1));

fvaltable(p,k,4) = out4.ferror;
errortable(p,k,4) = out4.error(end)/out1.error(1);

traintable(k,4,1) = norm(data.train_set1*out4.W(:,1) - data.train_label1)/norm(data.train_label1);
testtable(k,4,1) = norm(data.test_set1*out4.W(:,1) - data.test_label1)/norm(data.test_label1);
traintable(k,4,2) = norm(data.train_set2*out4.W(:,2) - data.train_label2)/norm(data.train_label2);
testtable(k,4,2) = norm(data.test_set2*out4.W(:,2) - data.test_label2)/norm(data.test_label2);
traintable(k,4,3) = norm(data.train_set3*out4.W(:,3) - data.train_label3)/norm(data.train_label3);
testtable(k,4,3) =  norm(data.test_set3*out4.W(:,3) - data.test_label3)/norm(data.test_label3);
traintable(k,4,4) = norm(data.train_set4*out4.W(:,4) - data.train_label4)/norm(data.train_label4);
testtable(k,4,4) = norm(data.test_set4*out4.W(:,4) - data.test_label4)/norm(data.test_label4);
traintable(k,4,5) = norm(data.train_set5*out4.W(:,5) - data.train_label5)/norm(data.train_label5);
testtable(k,4,5) = norm(data.test_set5*out4.W(:,5) - data.test_label5)/norm(data.test_label5);

figure(1)
semilogy([1 out1.error(1:opts.numall:end-10)/max(out2.error(1:opts.numall:end))],'r:.' ,'LineWidth',3)

hold on
semilogy( out2.error(1:opts.numall:end)/max(out2.error(1:opts.numall:end)),'b-','LineWidth',3)

semilogy(out3.error(1:opts.numall:end)/max(out2.error(1:opts.numall:end)),'k--','LineWidth',3) % Changed line style to '--' and LineWidth increased
semilogy(out4.error(1:opts.numall:end)/max(out2.error(1:opts.numall:end)),'m-.','LineWidth',3) % Changed line style to '-.' and LineWidth increased
%

figure(2)
semilogy((out1.fval(1:opts.numall:end))/abs(out1.fval(1)) ,'r:.' ,'LineWidth',3)
hold on
semilogy((out2.fval(1:opts.numall:end) )/abs(out1.fval(1)) ,'b-' ,'LineWidth',3)
semilogy((out3.fval(1:opts.numall:end) )/abs(out1.fval(1)) ,'k--' ,'LineWidth',3)
semilogy((out4.fval(1:opts.numall:end) )/abs(out1.fval(1)) ,'m-.' ,'LineWidth',3)



opts.options = 2;
opts.methods = 'sto';
opts.lr = @(i) min(max(6/i^(1/1.7),0.01),1);
opts.maxiter = 8000;
out1all = sto_ipm_multi_task3(blk,x0,data,At,b,opts);
out1all.errorend = min(out1all.error);
out1all.ferror = out1all.fval(end)/abs(out1.fval(1));

fvaltable(p,k,5) = out1all.ferror;
errortable(p,k,5) = out1all.error(end)/max(out1.error);


traintable(k,5,1) = norm(data.train_set1*out1all.W(:,1) - data.train_label1)/norm(data.train_label1);
testtable(k,5,1) = norm(data.test_set1*out1all.W(:,1) - data.test_label1)/norm(data.test_label1);
traintable(k,5,2) = norm(data.train_set2*out1all.W(:,2) - data.train_label2)/norm(data.train_label2);
testtable(k,5,2) = norm(data.test_set2*out1all.W(:,2) - data.test_label2)/norm(data.test_label2);
traintable(k,5,3) = norm(data.train_set3*out1all.W(:,3) - data.train_label3)/norm(data.train_label3);
testtable(k,5,3) =  norm(data.test_set3*out1all.W(:,3) - data.test_label3)/norm(data.test_label3);
traintable(k,5,4) = norm(data.train_set4*out1all.W(:,4) - data.train_label4)/norm(data.train_label4);
testtable(k,5,4) = norm(data.test_set4*out1all.W(:,4) - data.test_label4)/norm(data.test_label4);
traintable(k,5,5) = norm(data.train_set5*out1all.W(:,5) - data.train_label5)/norm(data.train_label5);
testtable(k,5,5) = norm(data.test_set5*out1all.W(:,5) - data.test_label5)/norm(data.test_label5);


opts.options = 3;
opts.methods = 'sto';
opts.lr = @(i) min(max(6/i^(1/1.7),0.01),1);
opts.maxiter = 4000;
out2all = sto_ipm_multi_task3(blk,x0,data,At,b,opts);
% out2all.errorend = min(out2all.error);
out2all.ferror = out2all.fval(end)/abs(out1.fval(1));
fvaltable(p,k,6) = out2all.ferror;
errortable(p,k,6) = out2all.error(end)/max(out1.error);

traintable(k,6,1) = norm(data.train_set1*out2all.W(:,1) - data.train_label1)/norm(data.train_label1);
testtable(k,6,1) = norm(data.test_set1*out2all.W(:,1) - data.test_label1)/norm(data.test_label1);
traintable(k,6,2) = norm(data.train_set2*out2all.W(:,2) - data.train_label2)/norm(data.train_label2);
testtable(k,6,2) = norm(data.test_set2*out2all.W(:,2) - data.test_label2)/norm(data.test_label2);
traintable(k,6,3) = norm(data.train_set3*out2all.W(:,3) - data.train_label3)/norm(data.train_label3);
testtable(k,6,3) =  norm(data.test_set3*out2all.W(:,3) - data.test_label3)/norm(data.test_label3);
traintable(k,6,4) = norm(data.train_set4*out2all.W(:,4) - data.train_label4)/norm(data.train_label4);
testtable(k,6,4) = norm(data.test_set4*out2all.W(:,4) - data.test_label4)/norm(data.test_label4);
traintable(k,6,5) = norm(data.train_set5*out2all.W(:,5) - data.train_label5)/norm(data.train_label5);
testtable(k,6,5) = norm(data.test_set5*out2all.W(:,5) - data.test_label5)/norm(data.test_label5);



figure(1)
set(gca, 'FontSize', 16) % Larger font size for axis ticks
xlabel('epochs', 'FontSize', 22) % Larger font size for x-axis label
ylabel('average relative stationary', 'FontSize', 22) % Larger font size for y-axis label
legend('SIPM-ME$^1$','SIPM-PM','SIPM-EM','SIPM-RM', 'FontSize', 18,  'Location','northeast','Interpreter','latex')  % Larrger font size for legend

save_path = ['./result/SDP_multitask_error', num2str(p),num2str(n0), num2str(opts.batchsize), '.png'];

saveas(gcf,save_path)


figure(2)
set(gca, 'FontSize', 16) % Larger font size for axis ticks
xlabel('epochs', 'FontSize', 22) % Larger font size for x-axis label
ylabel('relative objective', 'FontSize', 22) % Larger font size for y-axis label
legend('SIPM-ME$^1$','SIPM-PM','SIPM-EM','SIPM-RM', 'FontSize', 18,  'Location','northeast','Interpreter','latex')  % Lararger font size for legend

save_path = ['./result/SDP_multitask_ferror', num2str(p), num2str(n0),num2str(opts.batchsize), '.png'];

saveas(gcf,save_path)
filename = ['./result/SDP_multitask_out', num2str(p),num2str(n0),num2str(opts.batchsize),];
save(filename,"out1","out2","out3","out4",'out1all');
    end
end
%% nuclear norm
out5 = sto_nuclear_multi_task(blk,x0,data,At,b,opts);
out5.ferror = out5.fval(end) - opts.lambda2/2*trace(out5.W*out5.Omega{1}*out5.W');
out5.Omega{1} = out5.W'*out5.W/trace(out5.W'*out5.W);
out5.ferror = out5.ferror +  opts.lambda2/2*trace(out5.W*out5.Omega{1}*out5.W');

fvaltable(p,k,7) = out5.ferror/abs(out1.fval(1));
errortable(p,k,7) = out5.error(end)/max(out1.error);

traintable(k,7,1) = norm(data.train_set1*out5.W(:,1) - data.train_label1)/norm(data.train_label1);
testtable(k,7,1) = norm(data.test_set1*out5.W(:,1) - data.test_label1)/norm(data.test_label1);
traintable(k,7,2) = norm(data.train_set2*out5.W(:,2) - data.train_label2)/norm(data.train_label2);
testtable(k,7,2) = norm(data.test_set2*out5.W(:,2) - data.test_label2)/norm(data.test_label2);
traintable(k,7,3) = norm(data.train_set3*out5.W(:,3) - data.train_label3)/norm(data.train_label3);
testtable(k,7,3) =  norm(data.test_set3*out5.W(:,3) - data.test_label3)/norm(data.test_label3);
traintable(k,7,4) = norm(data.train_set4*out5.W(:,4) - data.train_label4)/norm(data.train_label4);
testtable(k,7,4) = norm(data.test_set4*out5.W(:,4) - data.test_label4)/norm(data.test_label4);
traintable(k,7,5) = norm(data.train_set5*out5.W(:,5) - data.train_label5)/norm(data.train_label5);
testtable(k,7,5) = norm(data.test_set5*out5.W(:,5) - data.test_label5)/norm(data.test_label5);

% figure(3)
figure('Position', [100, 100, 1600, 400]);  % 设置图形窗口的位置和大小

dataall = squeeze(traintable(2,:,:))'/2;
dataall = rearrangeColumns(dataall);
datasetstr = ["winequality-red","winequality-white","energy-efficency","airquality","abalone"];

for i = 1:5

    subplot('Position', [0.15+(i-1)*0.15 0.1 0.13 0.8]); 
    datatmp = dataall(i,:);
    hold on;  
    for j = 1:length(datatmp)
        bar(j, datatmp(j), 'DisplayName', ['Label ' num2str(j)]); 
    end
    hold off; 
    

    ylim([0, 0.5]);  
    xticks([]);  
    

    if i == 1
        ylabel('loss per task', 'FontSize', 18);  
        set(gca, 'YTick', 0:0.1:0.5);  
    else
        set(gca, 'YTick', []);  
    end
    

    title([datasetstr(i)], 'FontSize', 18);
end



legend('SIPM-ME$^1$','SIPM-ME$^+$','SIPM-PM','SIPM-EM','SIPM-RM','IPM-FG','AM', 'FontSize', 14,  'Location','best','Interpreter','latex')
saveas(gcf,'multitask_train.png')


figure(4)
dataall = squeeze(testtable(2,:,:))'/2;
dataall = rearrangeColumns(dataall);
datasetstr = ["winequality-red","winequality-white","energy-efficency","airquality","abalone"];
figure('Position', [100, 100, 1600, 400]); 


for i = 1:5

    subplot('Position', [0.15+(i-1)*0.15 0.1 0.13 0.8]); 
    datatmp = dataall(i,:);
    hold on;  
    for j = 1:length(datatmp)
        bar(j, datatmp(j), 'DisplayName', ['Label ' num2str(j)]); 
    
    end
    hold off; 

    ylim([0, 0.5]);  
    xticks([]); 
    

    if i == 1
        ylabel('loss per task', 'FontSize', 18); 
        set(gca, 'YTick', 0:0.1:0.5); 
    else
        set(gca, 'YTick', []);  
    end
    
    title([datasetstr(i)], 'FontSize', 18);
end



legend('SIPM-ME$^1$','SIPM-ME$^+$','SIPM-PM','SIPM-EM','SIPM-RM','IPM-FG','AM', 'FontSize', 14,  'Location','best','Interpreter','latex')
saveas(gcf,'multitask_test.png')

saveas(gcf,'multitask_test.png')


        filename = ['./result/SDP222_tabelout',char(probname), num2str(n0),num2str(opts.batchsize)];
        save(filename,"fvaltable","errortable","traintable","testtable");



function [AX, AXorg] = AXmap(X, K, At, Lchol)
AX = AXfun_sdpnal(K,At,X);
end

function Aty = Atymap(y, K, At, Lchol)
Aty = Atyfun_sdpnal(K, At, y);
end
