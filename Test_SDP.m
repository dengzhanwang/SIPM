clear
close all
probname = 'abalone';
probnameall = ["spambase", "covtype"];
n0all = [100 100; 500 50];
fvaltable = zeros(2,2,6);
errortable = zeros(2,2,6);
timetable = zeros(2,2,6);
% probname = 'spambase';
% probname = 'covtype';
% for test = 1:10
for k =1
    for j = 21
        % for [n0, batchsize] = n0all
        close all
        probname = probnameall(k);
        file = ['./rcp/' ,char(probname), '.data'];
        test_dir = './';
        addpath(genpath(test_dir));
        cluster = 10;
        n0 = n0all(j,1);
        opts.batchsize = n0all(j,2);
        batchnum = n0all(j,1)/n0all(j,2);
        [blk,At,C,b,W0] = rcpread_sto(file, cluster, n0);
        model.At = MatCell(At);
        model.b = b;
        model.C = MatCell(C);
        model.K = Cone.fromblk(blk);
        tmpC = C{1};
        rng(2023)

        Csto= zeros(size(W0,1),size(W0,1),1000);
        W1 = W0;
        for ii = 1:1000
            W1 = W1 + 0.01*randn(size(W1));
            Csto(:,:,ii) = -W1*W1';
        end
        opts.lr = 0.01;
        opts.scale_data = 1;
        trans = struct;

        opts.batch = 1:100:1000;

        opts.maxiter = 10;
        opts.kappa = 10;
        opts.gamma = 0.1;
        opts.lr = 1;


  
        At = model.At;
        K = model.K;
        b = model.b;
        C = model.C;
        opts.Amap = @(X) AXmap(X, blk, At);
        opts.ATmap = @(y) Atymap(y, blk, At);
        x0{1} = eye(n0)*(cluster)/n0 + (1-(cluster)/n0)/(n0-1)*(ones(n0) - eye(n0));
        %% solve
        opts.methods = 'sto';
        opts.options = 1;
        opts.lr = @(i) min(max(3/i^(1/2),0.01),1);
        out1 = sto_ipm_sdp(blk,x0,At,Csto,b,opts);
        out1.ferror = out1.fval(end)/abs(out1.fval(1));

        fvaltable(k,j,1) = out1.ferror;
        errortable(k,j,1) = out1.error(end)/out1.error(1);
        figure(1)
        semilogy(out1.error(1:batchnum:end)/out1.error(1),'r:.' ,'LineWidth',3) % LineWidth increased to 3 for thicker lines
        hold on
      
    end
end

filename = ['./result/SDP222_tabelout',char(probname), num2str(n0),num2str(opts.batchsize)];
save(filename,"fvaltable","errortable")

