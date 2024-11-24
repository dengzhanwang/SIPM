%% SOCP case1
clear
close all
seed = 2024;
rng(seed)
probname = 'covtype';

probnameall = ["winequality-white", "covtype"];

fname = [probname,'.data'];

test_dir = './';
addpath(genpath(test_dir));
n0all = [5000,500;10000,500];

fvaltable = zeros(10,2,7);
errortable = zeros(10,2,7);
timetable = zeros(10,2,7);

for test = 1
    for k = 1:2
        close all
        n0 = n0all(k,1);
        [blk,data,At] = Socp_energy_read(fname,n0);
        batchnum = n0all(k,1)/n0all(k,2);
        At{1}(11:end,:) = At{1}(11:end,:)/10;

        opts.batchsize = n0all(k,2);
        opts.lambda1 = 0.01;
        opts.lambda2 = 0.01;

        opts.maxiter = 200000;
        Amap = @(x) At'*x;
        ATmap = @(y) At*y;
        opts.Amap = Amap;
        opts.ATmap = ATmap;
        opts.n = size(At{1},2);
        opts.gamma = 0.1;
        x0 = zeros(opts.n*2+2,1);
        x0(1) = 10*test;
        x0(opts.n + 2) = 10*test;

        opts.options = 1;
        opts.methods = 'sto';
        opts.lr = @(i) min(max(3/i^(1/1.7),0.0001),1)*0.005;
        tt = tic;
        out1 =  sto_ipm_socp1(blk,x0,At,data,opts);
        out1.ferror = out1.fval(end)/abs(out1.fval(1));
        out1.time = toc(tt);

        fvaltable(test,k,1) = out1.ferror;
        errortable(test,k,1) = out1.error(end);
        timetable(test,k,1) = out1.time;

        figure(1)

        semilogy(out1.error(1:batchnum:end),'r:.' ,'LineWidth',3) % LineWidth increased to 3 for thicker lines
        out1.ferror = out1.fval(end)/abs(out1.fval(1));

       
    end
end

filename = ['./result/SOCP_table',char(probname), num2str(n0),num2str(opts.batchsize),];
save(filename,"fvaltable","errortable","timetable");
1;
%% alm-sto
% seed = 2024;
% rng(seed)
% [blk,data,At] = Socp_energy_read('ENB2012_data.xlsx');
% 
% 
% 
% %%
% clear
% clf
% [blk,data,At] = Socp_energy_read('ENB2012_data.xlsx');
% opts.n = size(At{1},2);
% opts.gamma = 0.1;
% opts.lambda1 = 0.01;
% opts.lambda2 = 0.01;
% opts.batchsize = 50;
% 
% data_train = data.train_set;
% label_train = data.train_label;
% data_test = data.test_set;
% label_test = data.test_label;
% x0 = zeros(opts.n*2+2,1);
% x0(1) = 1000;
% x0(opts.n + 2) = 1000;
% opts.maxiter = 10000;
% xnorm = [];
% reerror_train = [];
% reerror_test = [];
% for lambda = [0.0001 0.001 0.01 0.02 0.05 0.1 0.2]
%     opts.lambda1 = lambda;
%     opts.lambda2 = lambda;
%     opts.methods = 'mom';
%     opts.lr = @(i) min(max(6/i^(1/1.8),0.0001),1)*0.005;
%     out2 =  sto_ipm_socp1(blk,x0,At,data,opts);
%     xnorm = [xnorm norm(out2.x(2:9))];
%     reerror_train = [reerror_train  norm(data_train*out2.x(2:9) - label_train)/length(label_train)];
%     reerror_test = [reerror_test norm(data_test*out2.x(2:9) - label_test)/length(label_test)];
% end
% 
% data.train_miss = data.train_set;
% 
% xnorm_no = [];
% reerror_train_no = [];
% reerror_test_no = [];
% for lambda = [0.0001 0.001 0.01 0.02 0.05 0.1 0.2]
%     opts.lambda1 = lambda;
%     opts.lambda2 = lambda;
%     opts.methods = 'mom';
%     opts.lr = @(i) min(max(6/i^(1/1.8),0.0001),1)*0.005;
%     out3 =  sto_ipm_socp1(blk,x0,At,data,opts);
%     xnorm_no = [xnorm_no norm(out3.x(2:9))];
%     reerror_train_no = [reerror_train_no  norm(data_train*out3.x(2:9) - label_train)/length(label_train)];
%     reerror_test_no = [reerror_test_no norm(data_test*out3.x(2:9) - label_test)/length(label_test)];
% end
% 
% % figure(3)
% % dlambda = [0.001 0.01 0.02 0.05 0.1 0.2 0.5];
% %
% % % 绘制训练和测试误差
% % plot(dlambda, reerror_train, 'k-.', 'Marker', '*', 'LineWidth', 1.5)
% % hold on
% % plot(dlambda, reerror_test, 'k:.', 'Marker', '+', 'LineWidth', 1.5)
% %
% % plot(dlambda, reerror_train_no, 'b-.', 'Marker', '*', 'LineWidth', 1.5)
% % plot(dlambda, reerror_test_no, 'b:.', 'Marker', '+', 'LineWidth', 1.5)
% % ylim([0 0.5])
% % % 设置坐标轴标签和字体
% % xlabel('lambda ', 'FontWeight', 'bold', 'FontSize', 12)
% % ylabel('Regression error', 'FontWeight', 'bold', 'FontSize', 12)
% % legend('Train Error With Missing', 'Test Error With missing', 'Train Error ', 'Test Error ', 'Location', 'northwest')
% % % 设置坐标轴粗体
% % set(gca, 'FontWeight', 'bold', 'FontSize', 12)
% % saveas(gcf,'./result/SOCP_error1_comp.png')
% % hold off
% 
% 
% 1;
% % % 将表格数据转化为 double 类型矩阵
% % data_matrix = table2array(data);
% %% SOCP case2
% clear
% clf
% [blk,data,At] = Socp_energy_read('winequality-white.csv');
% % 读取红葡萄酒数据集
% 
% data_train = data.train_set;
% label_train = data.train_label;
% data_test = data.test_set;
% label_test = data.test_label;
% 
% opts.lambda1 = 0.01;
% opts.lambda2 = 0.01;
% opts.batchsize = 50;
% 
% opts.maxiter = 10000;
% Amap = @(x) At'*x;
% ATmap = @(y) At*y;
% opts.Amap = Amap;
% opts.ATmap = ATmap;
% opts.n = size(At{1},2);
% opts.gamma = 0.1;
% x0 = zeros(opts.n*2+2,1);
% x0(1) = 1000;
% x0(opts.n + 2) = 1000;
% 
% opts.methods = 'sto';
% opts.lr = @(i) min(max(6/i^(1/1.5),0.0001),1)*0.005;
% out1 =  sto_ipm_socp1(blk,x0,At,data,opts);
% figure(1)
% semilogy(out1.error,'r:.' ,'LineWidth',3) % LineWidth increased to 3 for thicker lines
% hold on
% 
% 
% opts.methods = 'mom';
% opts.lr = @(i) min(max(6/i^(1/1.8),0.0001),1)*0.005;
% out2 =  sto_ipm_socp1(blk,x0,At,data,opts);
% figure(1)
% semilogy(out2.error,'b-','LineWidth',3) % Changed line style to '-' and LineWidth increased
% 
% opts.gamma = 0.1;
% opts.methods = 'igt';
% opts.lr = @(i) min(max(6/i^(1/1.9),0.0001),1)*0.005;
% out3 =  sto_ipm_socp1(blk,x0,At,data,opts);
% figure(1)
% semilogy(out3.error,'k--','LineWidth',3) % Changed line style to '--' and LineWidth increased
% 
% opts.gamma = 0.1;
% opts.radius = 400000;
% opts.methods = 'recursiv';
% opts.lr = @(i) min(max(6/i^(1/2),0.0001),1)*0.005;
% out4 =  sto_ipm_socp1(blk,x0,At,data,opts);
% figure(1)
% semilogy(out4.error,'m-.','LineWidth',3) % Changed line style to '-.' and LineWidth increased
% 
% 
% figure(2)
% semilogy((out1.fval - min(out4.fval))/abs(out1.fval(1)) ,'r:.' ,'LineWidth',3)
% hold on
% semilogy((out2.fval - min(out4.fval))/abs(out1.fval(1)) ,'b-' ,'LineWidth',3)
% semilogy((out3.fval - min(out4.fval))/abs(out1.fval(1)) ,'k--' ,'LineWidth',3)
% semilogy((out4.fval - min(out4.fval)+0.1)/abs(out1.fval(1)) ,'m-.' ,'LineWidth',3)
% 
% % Set font size for labels and legend
% figure(1)
% set(gca, 'FontSize', 16) % Larger font size for axis ticks
% xlabel('Iteration', 'FontSize', 18) % Larger font size for x-axis label
% ylabel('Error', 'FontSize', 18) % Larger font size for y-axis label
% legend('SIPM-ME','SIPM-PM','SIPM-EM','SIPM-RM', 'FontSize', 14, 'Location','northeast') % Larger font size for legend
% saveas(gcf,'./result/SOCP2_error.png')
% 
% figure(2)
% set(gca, 'FontSize', 16) % Larger font size for axis ticks
% xlabel('Iteration', 'FontSize', 18) % Larger font size for x-axis label
% ylabel('$(f-f^*)/f^*$', 'FontSize', 18,'Interpreter','latex') % Larger font size for y-axis label
% legend('SIPM-ME','SIPM-PM','SIPM-EM','SIPM-RM', 'FontSize', 14, 'Location','northeast') % Larger font size for legend
% saveas(gcf,'./result/SOCP2_ferror.png')
% 1;
% 
% % test
% xbest = pinv(data_train)*label_train;
% norm(data_train*xbest - label_train)/length(label_train)
% norm(xbest)
% norm(data_train*out1.x(2:9) - label_train)/length(label_train)
% norm(out4.x(2:9))
% 
% norm(data_test*xbest - label_test)^2/length(label_test)
% norm(data_test*out1.x(2:9) - label_test)^2/length(label_test)
% %%
% clear
% clf
% [blk,data,At] = Socp_energy_read('winequality-red.csv');
% opts.n = size(At{1},2);
% opts.gamma = 0.1;
% opts.lambda1 = 0.01;
% opts.lambda2 = 0.01;
% opts.batchsize = 50;
% 
% data_train = data.train_set;
% label_train = data.train_label;
% data_test = data.test_set;
% label_test = data.test_label;
% x0 = zeros(opts.n*2+2,1);
% x0(1) = 1000;
% x0(opts.n + 2) = 1000;
% opts.maxiter = 10000;
% xnorm = [];
% reerror_train = [];
% reerror_test = [];
% for lambda = [0.001 0.01 0.02 0.05 0.1 0.2 0.5]
%     opts.lambda1 = 0.001;
%     opts.lambda2 = lambda;
%     opts.methods = 'mom';
%     opts.lr = @(i) min(max(6/i^(1/1.8),0.0001),1)*0.005;
%     out2 =  sto_ipm_socp1(blk,x0,At,data,opts);
%     xnorm = [xnorm norm(out2.x(2:9))];
%     reerror_train = [reerror_train  norm(data_train*out2.x(2:9) - label_train)/length(label_train)];
%     reerror_test = [reerror_test norm(data_test*out2.x(2:9) - label_test)/length(label_test)];
% end
% 
% data.train_miss = data.train_set;
% 
% xnorm_no = [];
% reerror_train_no = [];
% reerror_test_no = [];
% for lambda = [0.001 0.01 0.02 0.05 0.1 0.2 0.5]
%     opts.lambda1 = 0.001;
%     opts.lambda2 = lambda;
%     opts.methods = 'mom';
%     opts.lr = @(i) min(max(6/i^(1/1.8),0.0001),1)*0.005;
%     out3 =  sto_ipm_socp1(blk,x0,At,data,opts);
%     xnorm_no = [xnorm_no norm(out3.x(2:9))];
%     reerror_train_no = [reerror_train_no  norm(data_train*out3.x(2:9) - label_train)/length(label_train)];
%     reerror_test_no = [reerror_test_no norm(data_test*out3.x(2:9) - label_test)/length(label_test)];
% end
% 
% figure(3)
% dlambda = [0.001 0.01 0.02 0.05 0.1 0.2 0.5];
% 
% % 绘制训练和测试误差
% plot(dlambda, reerror_train, 'k-.', 'Marker', '*', 'LineWidth', 1.5)
% hold on
% plot(dlambda, reerror_test, 'k:.', 'Marker', '+', 'LineWidth', 1.5)
% 
% plot(dlambda, reerror_train_no, 'b-.', 'Marker', '*', 'LineWidth', 1.5)
% plot(dlambda, reerror_test_no, 'b:.', 'Marker', '+', 'LineWidth', 1.5)
% 
% % 设置坐标轴标签和字体
% xlabel('lambda ', 'FontWeight', 'bold', 'FontSize', 12)
% ylabel('Regression error', 'FontWeight', 'bold', 'FontSize', 12)
% legend('Train Error With Missing', 'Test Error With missing', 'Train Error ', 'Test Error ', 'Location', 'southeast')
% % 设置坐标轴粗体
% set(gca, 'FontWeight', 'bold', 'FontSize', 12)
% saveas(gcf,'./result/SOCP_error2_comp.png')
% hold off