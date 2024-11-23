
clear
close all
seed = 2024;
rng(seed)
probname = 'winequality-white';


test_dir = './';
addpath(genpath(test_dir));
n0all = [2000,200;4000,200];

   fname = [probname,'.csv'];

fvaltable = zeros(10,2,6);
errortable = zeros(10,2,6);
timetable = zeros(10,2,6);
for test = 1
    for k = 1:2
        close all
        n0 = n0all(k,1);
        [blk,data,At] = Socp_energy_read(fname,n0);
        At{1}(11:end,:) = At{1}(11:end,:)/100;
        opts.batchsize = n0all(k,2);
        opts.lambda1 = 0.01;
        opts.lambda2 = 0.01;
        batchnum = n0all(k,1)/n0all(k,2);
        opts.maxiter = 200000;
        Amap = @(x) At'*x;
        ATmap = @(y) At*y;
        opts.Amap = Amap;
        opts.ATmap = ATmap;
        opts.n = size(At{1},2);
        opts.gamma = 0.1;
        x0 = zeros(opts.n*2+2,1);
        x0(1) = 1000*test;
        x0(opts.n + 2) = 1000*test;

        opts.options = 1;
        opts.methods = 'sto';
        opts.lr = @(i) min(max(3/i^(1/1.6),0.0001),1)*0.005;
        tt = tic;
        out1 =  sto_ipm_socp1(blk,x0,At,data,opts);
        out1.ferror = out1.fval(end)/abs(out1.fval(1));
        out1.time = toc(tt);

        fvaltable(test,k,1) = out1.ferror;
        errortable(test,k,1) = out1.error(end);
        timetable(test,k,1) = out1.time;

        figure(1)

        semilogy(out1.error(1:batchnum:end),'r:.' ,'LineWidth',3) 
        out1.ferror = out1.fval(end)/abs(out1.fval(1));

        hold on


        opts.methods = 'mom';
        opts.lr = @(i) min(max(4/i^(1/1.7),0.0001),1)*0.005;
        tt = tic;
        out2 =  sto_ipm_socp1(blk,x0,At,data,opts);
        out2.ferror = out2.fval(end)/abs(out1.fval(1));
        out2.time = toc(tt);
        fvaltable(test,k,2) = out2.ferror;
        errortable(test,k,2) = out2.error(end);
        timetable(test,k,2) = out2.time;
        figure(1)
        semilogy(out2.error(1:batchnum:end),'b-','LineWidth',3) 
        opts.gamma = 0.1;
        opts.methods = 'igt';
        opts.lr = @(i) min(max(4/i^(1/1.6),0.0001),1)*0.005;
        tt = tic;
        out3 =  sto_ipm_socp1(blk,x0,At,data,opts);
        out3.ferror = out3.fval(end)/abs(out1.fval(1));
        out3.time = toc(tt);

        fvaltable(test,k,3) = out3.ferror;
        errortable(test,k,3) = out3.error(end);
        timetable(test,k,3) = out3.time;
        figure(1)
        semilogy(out3.error(1:batchnum:end),'k--','LineWidth',3) 
        %
        opts.gamma = 0.1;
        opts.radius = 4e10;
        opts.methods = 'recursiv';
        opts.lr = @(i) min(max(7/i^(1/1.7),0.0001),1)*0.005;
        tt = tic;
        out4 =  sto_ipm_socp1(blk,x0,At,data,opts);
out4.time = toc(tt);
        out4.ferror = out4.fval(end)/abs(out1.fval(1));

        fvaltable(test,k,4) = out4.ferror;
        errortable(test,k,4) = out4.error(end);
        timetable(test,k,4) = out4.time;
        
        figure(1)
        semilogy(out4.error(1:batchnum:end),'m-.','LineWidth',3) % Changed line style to '-.' and LineWidth increased

        opts.options = 2;

        opts.methods = 'sto';
        opts.lr = @(i) min(max(4/i^(1/1.7),0.0001),1)*0.005;
        tt = tic;
        out1all =  sto_ipm_socp1(blk,x0,At,data,opts);
        out1all.time = toc(tt);
        out1all.ferror = out1all.fval(end)/abs(out1all.fval(1));
        fvaltable(test,k,5) = out1all.ferror;
        errortable(test,k,5) = out1all.error(end);
        timetable(test,k,5) = out1all.time;

        figure(2)
        semilogy((out1.fval(1:batchnum:end))/abs(out1.fval(1)) ,'r:.' ,'LineWidth',3)
        hold on
        semilogy((out2.fval(1:batchnum:end) )/abs(out1.fval(1)) ,'b-' ,'LineWidth',3)
        semilogy((out3.fval(1:batchnum:end))/abs(out1.fval(1)) ,'k--' ,'LineWidth',3)
        semilogy((out4.fval(1:batchnum:end))/abs(out1.fval(1)) ,'m-.' ,'LineWidth',3)

        % Set font size for labels and legend
        if opts.batchsize == 200 && n0 == 4000
            figure(1)
            set(gca, 'FontSize', 16) % Larger font size for axis ticks
            xlabel('epochs', 'FontSize', 22) % Larger font size for x-axis label
            ylabel('average relative stationary', 'FontSize', 22) % Larger font size for y-axis label
            legend('SIPM-ME$^1$','SIPM-PM','SIPM-EM','SIPM-RM', 'FontSize', 18,  'Location','northeast','Interpreter','latex')  % Larger font size for legend
            save_path = ['./result/SOCP_error', num2str(n0),char(probname), num2str(opts.batchsize), '.png'];
            saveas(gcf,save_path)

            figure(2)
            set(gca, 'FontSize', 16) % Larger font size for axis ticks
            xlabel('epochs', 'FontSize', 22) % Larger font size for x-axis label
            ylabel('relative objective', 'FontSize', 22) % Larger font size for y-axis label
            legend('SIPM-ME$^1$','SIPM-PM','SIPM-EM','SIPM-RM', 'FontSize', 18,  'Location','northeast','Interpreter','latex')  % Larger font size for legend
            % legend('sto','mom','igt','recursiv', 'FontSize', 14, 'Location','northeast') % Larger font size for legend
            save_path = ['./result/SOCP_ferror', num2str(n0),char(probname), num2str(opts.batchsize), '.png'];
            saveas(gcf,save_path)
        end

        opts2.rho = 1;
        opts2.eta = 0.01;
        opts2.alpha = 0.9;
        opts2.n = 10;
        opts2.theta = 1;
        opts.n = size(At{1},2);

        opts2.lambda1 = 0.01;
        opts2.lambda2 = 0.01;
        opts2.batchsize = 500;
        opts2.max_iter = 10000;
        opts2.n = size(At{1},2);
        opts2.gamma = 0.1;

        tt = tic;
        out_al =  alm_socp1(x0,At,data,opts2);
        out_al.time = toc(tt);
        out_al.ferror = out_al.fval(end)/abs(out1.fval(1));
        fvaltable(test,k,6) = out_al.ferror;
        timetable(test,k,6) = out_al.time;

        filename = ['./result/SOCP_out',char(probname), num2str(n0),num2str(opts.batchsize),];
        save(filename,"out1","out2","out3","out4",'out1all',"out_al");
    end
end

filename = ['./result/SOCP_table',char(probname), num2str(n0),num2str(opts.batchsize),];
save(filename,"fvaltable","errortable","timetable");
1;