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
    for j = 2
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
        % if n0 == 500
        [U,S,V] = svd(max(out1.x{1},0));
        dS = sqrt(diag(S));
        figure(3)
        rankx = 5
        U2 = U(:,1:rankx);
        Wnew = U2.*dS(1:rankx)';
        idx = kmeans(Wnew,3);
        idx2 = kmeans(W0,3);
        [U1,S1,V1] = svd(W0);
        dS1 = U1*S1(:,1:3)/1000;
        for ii = 1:1000
            dS1 = dS1 + 0.005*randn(size(dS1));
            % Csto(:,:,ii) = -W1*W1';

            colors = [1 0 0; 1 0 1;0 0 1];


            % 根据 idx 中的值从 colors 中选择颜色
            scatterColors = colors(idx, :);
            figure(3)
            scatter3(dS1(:,1), dS1(:,2), dS1(:,3),50, scatterColors, 'filled')
            % xlim([-10 0])
            % ylim([-4 4])
            % zlim([-1 2])
            xlim([-4 1])
            ylim([-2 0])
            zlim([-0.4 0.6])
            set(gca, 'FontSize', 16) % Larger font size for axis ticks
            if strcmp(char(probname),'spambase')
                view(10,10)
            else
                view(80,-10)
            end

            if ii == 1 || ii== 333 || ii == 666 || ii == 1000
                figure(3)
                save_path = ['./result/SDP222_cluster',num2str(ii), num2str(n0),char(probname),num2str(opts.batchsize), '.png'];
                saveas(gcf,save_path)
            end
        end
        matchedIndex = matchClusterLabels(idx, idx2);
        figure(4)
        scatterColors2 = colors(matchedIndex , :);
        scatter3(dS1(:,1), dS1(:,2), dS1(:,3),50, scatterColors2, 'filled')

        set(gca, 'FontSize', 16) % Larger font size for axis ticks

        colormap(jet)
        hold on



        figure(4)
        save_path = ['./result/SDP222_cluster_ori', num2str(n0),char(probname),num2str(opts.batchsize), '.png'];
        saveas(gcf,save_path)



        opts.methods = 'mom';
        opts.lr = @(i) min(max(3/i^(1/2),0.01),12);
        out2 = sto_ipm_sdp(blk,x0,At,Csto,b,opts);
        out2.ferror = out2.fval(end)/abs(out1.fval(1));

        fvaltable(k,j,2) = out2.ferror;
        errortable(k,j,2) = out2.error(end);
        out2.errorend = min(out2.error)/out1.error(1);

        figure(1)
        semilogy(out2.error(1:batchnum:end)/out1.error(1),'b-','LineWidth',3) % Changed line style to '-' and LineWidth increased



        opts.methods = 'igt';
        opts.lr = @(i) min(max(4/i^(1/2),0.01),1);
        out3 = sto_ipm_sdp(blk,x0,At,Csto,b,opts);
        out3.ferror = out3.fval(end)/abs(out1.fval(1));

        fvaltable(k,j,3) = out3.ferror;
        errortable(k,j,3) = out3.error(end);
        out3.errorend = min(out3.error);
        figure(1)
        semilogy(out3.error(1:batchnum:end)/out1.error(1),'k--','LineWidth',3) % Changed line style to '--' and LineWidth increased


        hold on

        opts.radius = 10e12;
        opts.methods = 'recursiv';
        opts.lr = @(i) min(max(5/i^(1/2),0.01),1);
        out4 = sto_ipm_sdp(blk,x0,At,Csto,b,opts);

        out4.ferror = out4.fval(end)/abs(out1.fval(1));

        fvaltable(k,j,4) = out4.ferror;
        errortable(k,j,4) = out4.error(end);
        out4.errorend = min(out4.error);
        figure(1)
        semilogy(out4.error(1:batchnum:end)/out1.error(1),'m-.','LineWidth',3) % Changed line style to '-.' and LineWidth increased

        opts.methods = 'sto';
        opts.options = 2;
        opts.lr = @(i) min(max(3/i^(1/2),0.01),1);
        out1all = sto_ipm_sdp(blk,x0,At,Csto,b,opts);
        out1all.errorend = min(out1all.error);
        out1all.ferror = out1all.fval(end)/abs(out1.fval(1));

        fvaltable(k,j,5) = out1all.ferror;
        errortable(k,j,5) = out1all.error(end);

        opts.methods = 'sto';
        opts.options = 3;
        opts.lr = @(i) min(max(3/i^(1/2),0.01),1);
        out2all = sto_ipm_sdp(blk,x0,At,Csto,b,opts);
        out2all.errorend = min(out1all.error);
        out2all.ferror = out1all.fval(end)/abs(out1.fval(1));

        fvaltable(k,j,6) = out2all.ferror;
        errortable(k,j,6) = out2all.error(end);


        figure(2)
        plot((out1.fval(1:batchnum:end))/abs(out1.fval(1) )+2 ,'r:.' ,'LineWidth',3)
        hold on
        plot(out2.fval(1:batchnum:end)/abs(out1.fval(1))+2 ,'b-' ,'LineWidth',3)
        plot(out3.fval(1:batchnum:end)/abs(out1.fval(1))+2 ,'k--' ,'LineWidth',3)
        plot(out4.fval(1:batchnum:end)/abs(out1.fval(1))+2 ,'m-.' ,'LineWidth',3)
        low = 1.2*min(out4.fval(1:out4.outiter)/abs(out1.fval(1)));
        ylim([low -1])
        % if n0 == 500
        % Set font size for labels and legend
        figure(1)
        set(gca, 'FontSize', 16) % Larger font size for axis ticks
        xlabel('epochs', 'FontSize', 22) % Larger font size for x-axis label
        ylabel('average relative stationary', 'FontSize', 22) % Larger font size for y-axis label
        legend('SIPM-ME$^1$','SIPM-PM','SIPM-EM','SIPM-RM', 'FontSize', 18,  'Location','northeast','Interpreter','latex')  % Larger font size for legend

        save_path = ['./result/SDP222_error', num2str(n0),char(probname), num2str(opts.batchsize), '.png'];

        saveas(gcf,save_path)


        figure(2)
        set(gca, 'FontSize', 16) % Larger font size for axis ticks
        xlabel('epochs', 'FontSize', 22) % Larger font size for x-axis label
        ylabel('relative objective', 'FontSize', 22) % Larger font size for y-axis label
        legend('SIPM-ME$^1$','SIPM-PM','SIPM-EM','SIPM-RM', 'FontSize', 18,  'Location','northeast','Interpreter','latex')  % Larger font size for legend

        save_path = ['./result/SDP222_ferror', num2str(n0),char(probname),num2str(opts.batchsize), '.png'];

        saveas(gcf,save_path)
        filename = ['./result/SDP_out',char(probname), num2str(n0),num2str(opts.batchsize),];
        save(filename,"out1","out2","out3","out4",'out1all');
        1;

        1;
    end
end

filename = ['./result/SDP222_tabelout',char(probname), num2str(n0),num2str(opts.batchsize)];
save(filename,"fvaltable","errortable")
% end
%% LP case Sigmod
% clear
% clf
% rng(2025)
% [blk, At,Q, b ] = Lp_liver_Sigmod_read('liver/bupa');
% opts.Amap = @(X) AXmap(X, blk, At);
% opts.ATmap = @(y) Atymap(y, blk, At);
%
% trans = struct;
%
% opts.batch = 1:100:1000;
% opts.batchsize = 100;
% opts.maxiter = 1000;
% opts.kappa = 0.01;
% opts.gamma = 0.1;
% opts.radius = 40;
%
% opts.lr = 1;
% pindex = find(At{1}(:,1) == 1);
% nindex = find(At{1}(:,1) == -1);
% n0 = length(b) - 1;
% x(pindex) = length(nindex) / length(pindex);
% x(nindex) = 1;
% x(n0+1:2*n0) = 5 -x(1:n0);
% x0{1} = x';
% % solve
% opts.methods = 'sto';
% opts.lr = @(i) min(max(4/i^(3/4),0.01),1);
% out1 = sto_ipm_sdp(blk,x0,At,Q,b,opts);
% figure(1)
% semilogy(out1.error,'r:.' ,'LineWidth',3) % LineWidth increased to 3 for thicker lines
% hold on
%
% % opts.maxiter = 1000;
% opts.methods = 'mom';
% opts.lr = @(i) min(max(8/i^(3/4),0.01),1);
% out2 = sto_ipm_sdp(blk,x0,At,Q,b,opts);
%
% figure(1)
% semilogy(out2.error,'b-','LineWidth',3) % Changed line style to '-' and LineWidth increased
%
%
% opts.gamma = 0.1;
% opts.methods = 'igt';
% opts.lr = @(i) min(max(8/i^(3/4),0.01),1);
% out3 = sto_ipm_sdp(blk,x0,At,Q,b,opts);
%
% figure(1)
% semilogy(out3.error,'k--','LineWidth',3) % Changed line style to '--' and LineWidth increased
%
%
% hold on
% opts.gamma = 0.1;
% opts.radius = 400000;
% opts.methods = 'recursiv';
% opts.lr = @(i) min(max(6/i^(2.8/4),0.01),1);
% out4 = sto_ipm_sdp(blk,x0,At,Q,b,opts);
% figure(1)
% semilogy(out4.error,'m-.','LineWidth',3) % Changed line style to '-.' and LineWidth increased
%
%
% figure(2)
% semilogy((out1.fval - min(out4.fval))/abs(out1.fval(1)) ,'r:.' ,'LineWidth',3)
% hold on
% semilogy((out2.fval - min(out4.fval))/abs(out1.fval(1)) ,'b-' ,'LineWidth',3)
% semilogy((out3.fval - min(out4.fval))/abs(out1.fval(1)) ,'k--' ,'LineWidth',3)
% semilogy((out4.fval - min(out4.fval))/abs(out1.fval(1)) ,'m-.' ,'LineWidth',3)
%
% % Set font size for labels and legend
% figure(1)
% set(gca, 'FontSize', 16) % Larger font size for axis ticks
% xlabel('Iteration', 'FontSize', 18) % Larger font size for x-axis label
% ylabel('Error', 'FontSize', 18) % Larger font size for y-axis label
% legend('sto','mom','igt','recursiv', 'FontSize', 14) % Larger font size for legend
% saveas(gcf,'./result/LP_Sigmod_error.png')
%
% figure(2)
% set(gca, 'FontSize', 16) % Larger font size for axis ticks
% xlabel('Iteration', 'FontSize', 18) % Larger font size for x-axis label
% ylabel('$(f-f^*)/f^*$', 'FontSize', 18,'Interpreter','latex') % Larger font size for y-axis label
% legend('sto','mom','igt','recursiv', 'FontSize', 14, 'Location','best') % Larger font size for legend
% saveas(gcf,'./result/LP_Sigmod_ferror.png')
% 1;
% %% LP case Gaussian
% clear
% clf
% rng(2025)
% [blk, At,Q, b ] = Lp_liver_Gaussian_read('liver/bupa');
% opts.Amap = @(X) AXmap(X, blk, At);
% opts.ATmap = @(y) Atymap(y, blk, At);
%
% trans = struct;
%
% opts.batch = 1:100:1000;
% opts.batchsize = 100;
% opts.maxiter = 1000;
% opts.kappa = 0.01;
% opts.gamma = 0.1;
% opts.radius = 40;
%
% opts.lr = 1;
% pindex = find(At{1}(:,1) == 1);
% nindex = find(At{1}(:,1) == -1);
% n0 = length(b) - 1;
% x(pindex) = length(nindex) / length(pindex);
% x(nindex) = 1;
% x(n0+1:2*n0) = 5 -x(1:n0);
% x0{1} = x';
% % solve
% opts.methods = 'sto';
% opts.lr = @(i) min(max(2.5/i^(3/4),0.01),1);
% out1 = sto_ipm_sdp(blk,x0,At,Q,b,opts);
% figure(1)
% semilogy(out1.error,'r:.' ,'LineWidth',3) % LineWidth increased to 3 for thicker lines
% hold on
%
%
% % opts.maxiter = 1000;
% opts.methods = 'mom';
% opts.lr = @(i) min(max(8/i^(3/4),0.01),1);
% out2 = sto_ipm_sdp(blk,x0,At,Q,b,opts);
%
% figure(1)
% semilogy(out2.error,'b-','LineWidth',3) % Changed line style to '-' and LineWidth increased
%
% opts.gamma = 0.1;
% opts.methods = 'igt';
% opts.lr = @(i) min(max(4/i^(3/4),0.01),1);
% out3 = sto_ipm_sdp(blk,x0,At,Q,b,opts);
%
% figure(1)
% semilogy(out3.error,'k--','LineWidth',3) % Changed line style to '--' and LineWidth increased
%
%
%
%
% hold on
% opts.gamma = 0.5;
% opts.radius = 400000;
% opts.methods = 'recursiv';
% opts.lr = @(i) min(max(4/i^(3/4),0.01),1);
% out4 = sto_ipm_sdp(blk,x0,At,Q,b,opts);
% figure(1)
% semilogy(out4.error,'m-.','LineWidth',3) % Changed line style to '-.' and LineWidth increased
%
%
% figure(2)
% semilogy((out1.fval - min(out4.fval))/abs(out1.fval(1)) ,'r:.' ,'LineWidth',3)
% hold on
% semilogy((out2.fval - min(out4.fval))/abs(out1.fval(1)) ,'b-' ,'LineWidth',3)
% semilogy((out3.fval - min(out4.fval))/abs(out1.fval(1)) ,'k--' ,'LineWidth',3)
% semilogy((out4.fval - min(out4.fval))/abs(out1.fval(1)) ,'m-.' ,'LineWidth',3)
%
% % Set font size for labels and legend
% figure(1)
% set(gca, 'FontSize', 16) % Larger font size for axis ticks
% xlabel('Iteration', 'FontSize', 18) % Larger font size for x-axis label
% ylabel('Error', 'FontSize', 18) % Larger font size for y-axis label
% legend('sto','mom','igt','recursiv', 'FontSize', 14) % Larger font size for legend
% saveas(gcf,'./result/LP_Gaussian_error.png')
%
% figure(2)
% set(gca, 'FontSize', 16) % Larger font size for axis ticks
% xlabel('Iteration', 'FontSize', 18) % Larger font size for x-axis label
% ylabel('$(f-f^*)/f^*$', 'FontSize', 18,'Interpreter','latex') % Larger font size for y-axis label
% legend('sto','mom','igt','recursiv', 'FontSize', 14, 'Location','northeast') % Larger font size for legend
% saveas(gcf,'./result/LP_Gaussian_ferror.png')
% 1;




% function [model_new,trans] = preprocess_SDP(model,opts,trans)
% K = model.K;
% At = model.At;
% b = model.b;
% C = model.C;
%
% %% scale the data and cholesky decomposition
% trans.mdim = length(b);
% trans.nblock = length(K);
% DA = speye(trans.mdim);
% bscale = 1;
% Cscale = 1;
%
% if opts.scale_data
%     % normA: norm of each row of A
%     normA = zeros(trans.mdim,1);
%     for k = 1:(trans.nblock)
%         normA = normA + sum(At{k}.*At{k})';
%     end
%     normA = max(1,sqrt(normA));
%     DA = spdiag(1./normA);
% end
% for k = 1:trans.nblock
%     At{k} = At{k}*DA;
% end
% b = DA*b;
% scale.DA = DA;
%
% AAt = AAtfun(At);
% m = size(AAt,1);
% Lchol = struct;
% % if trans.use_AAtchol
% %     if (nnz(AAt) < 0.2*m*m); use_chol=1; else,  use_chol=0; end
% % else
% use_chol= 0;
% % end
% if (use_chol)
%     [Lchol.R, Lchol.p, Lchol.perm] = chol(AAt,'vector');
%     Lchol.Rt = Lchol.R';
%     Lchol.matfct_options = 'spcholmatlab';
% else
%     if issparse(AAt); AAt = full(AAt); end;
%     Lchol.matfct_options = 'chol';
%     Lchol.perm = 1:m;
%     %     [Lchol.R, indef] = chol(AAt);
%     Lchol.R = speye(size(AAt));
%     Lchol.Rt = Lchol.R';
%     Lchol.isidentity = true;
% end
%
% % if isequal(Lchol.R ,eye(size(Lchol.R)),'fro') < 1e-10
% if isdiag(Lchol.R) && (norm(diag(Lchol.R) - 1) < 1e-7)
%     Lchol.isidentity = true;
% else
%     Lchol.isidentity = false;
% end
% b = fwsolve(Lchol,b);
%
% if (opts.scale_data==1)
%     bscale = max(1,norm(b));
%     Cscale = max(1,norm(C));
% end
%
% b = b / bscale; %单位化
% C = 1 / Cscale * C;
%
% model_new.K = K;
% model_new.At = At;
% model_new.b = b;
% model_new.C = C;
%
% objscale = bscale*Cscale;
% scale.bscale = bscale;
% scale.Cscale = Cscale;
% scale.objscale = objscale;
% scale.scale_data = opts.scale_data;
% trans.scale = scale;
% trans.Lchol = Lchol;
% end
%
function [AX, AXorg] = AXmap(X, K, At, Lchol)
AX = AXfun_sdpnal(K,At,X);
% AX = fwsolve(Lchol, AXorg);
end

function Aty = Atymap(y, K, At, Lchol)
Aty = Atyfun_sdpnal(K, At, y);
% Aty = Atyfun(K, At, bwsolve(Lchol, y));
end
