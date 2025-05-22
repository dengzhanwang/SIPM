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
