function out = sto_ipm_multi_task(blk,x0,data,At,b,opts)
par.numcolAt = length(b);
% kappa = opts.kappa;
% C{1} = Csto(:,:,1);
W= 0.01*ones(size(x0));
Omega{1} = opts.Omega0;
z = Omega;

[At,~,Omega,z,par.permA,par.permZ] = sortA(blk,At,[],b,Omega,z);
par.spdensity   = 0.4;
par.smallblkdim = 50;
[par.isspA,par.nzlistA,par.nzlistAsum,par.isspAy,par.nzlistAy] = nzlist(blk,At,par);
batchsize = opts.batchsize;
muk = 1;
trans.Amap = opts.Amap;
trans.ATmap = opts.ATmap;
gamma = opts.gamma;
lr = opts.lr;
x = x0;
Omegaold = Omega;
Wold = W;
out.error = [];
out.fval = [];
numall = size(data.train_set1,1)/batchsize;
% gradientave = sum(Csto(:,:,:),3)/batchsize;
for i = 1:opts.maxiter
    if mod(i,5) + 1 == 2
        train_set = data.train_set1;
        train_label = data.train_label1;
        batchsize = length(train_label);
    elseif mod(i,5) + 1 == 3
        train_set = data.train_set2;
        train_label = data.train_label2;
        batchsize = length(train_label);
    elseif mod(i,5) + 1 == 4
        train_set = data.train_set3;
        train_label = data.train_label3;
        batchsize = length(train_label);
    elseif mod(i,5) + 1 == 5
        train_set = data.train_set4;
        train_label = data.train_label4;
        batchsize = length(train_label);
    elseif mod(i,5) + 1 == 1
        train_set = data.train_set5;
        train_label = data.train_label5;
        batchsize = length(train_label);
    end
    % batchnum = mod(i,10);
    % batchnum = 0;
    % batch = batchnum*batchsize + 1: batchsize*(batchnum+1);
    if strcmp(opts.methods,'sto')
        [U, S] = eig(Omega{1});
        invomega = U * diag(1./(diag(S)+0.01)) * U';
        tmpW = W*invomega;
        % invomega = inv(Omega);
        if mod(i,5) + 1 == 2
            tmpgradient = train_set'*(train_set*W(:,1)-train_label); 
            gradient1 = 2*tmpgradient./(1+tmpgradient.^2)/batchsize*1e4;
            gradient1 = gradient1 + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1);
            gradient1 = [gradient1, zeros(size(gradient1)) ,zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 3
            % gradient1 = train_set'*(train_set*W(:,2)-train_label)/batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2);
                   tmpgradient = train_set'*(train_set*W(:,2)-train_label); 
            gradient1 = 2*tmpgradient./(1+tmpgradient.^2)/batchsize*1e4;
            gradient1 = gradient1 + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2);
            gradient1 = [zeros(size(gradient1)), gradient1,  zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 4
            % gradient1 = train_set'*(train_set*W(:,3)-train_label)/batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3);
                          tmpgradient = train_set'*(train_set*W(:,3)-train_label); 
            gradient1 = 2*tmpgradient./(1+tmpgradient.^2)/batchsize*1e4;
            gradient1 = gradient1 + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), gradient1,zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 5
            % gradient1 = train_set'*(train_set*W(:,4)-train_label)/batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4);
                          tmpgradient = train_set'*(train_set*W(:,4)-train_label); 
            gradient1 = 2*tmpgradient./(1+tmpgradient.^2)/batchsize*1e4;
            gradient1 = gradient1 + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1)), gradient1, zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 1
            % gradient1 = train_set'*(train_set*W(:,5)-train_label)/batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5);
                          tmpgradient = train_set'*(train_set*W(:,5)-train_label); 
            gradient1 = 2*tmpgradient./(1+tmpgradient.^2)/batchsize*1e4;
            gradient1 = gradient1 + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), zeros(size(gradient1)),zeros(size(gradient1)),gradient1];
        end
        gradient =  -opts.lambda2* invomega*W'*W*invomega   - muk * invomega;
    elseif strcmp(opts.methods,'mom')
        [U, S] = eig(Omega{1});
        invomega = U * diag(1./(diag(S)+0.01)) * U';
        tmpW = W*invomega;
        % invomega = inv(Omega);
        if mod(i,5) + 1 == 2
            gradient1 = train_set'*(train_set*W(:,1)-train_label)/batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1);
            gradient1 = [gradient1, zeros(size(gradient1)) ,zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 3
            gradient1 = train_set'*(train_set*W(:,2)-train_label)/batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2);
            gradient1 = [zeros(size(gradient1)), gradient1,  zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 4
            gradient1 = train_set'*(train_set*W(:,3)-train_label)/batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), gradient1,zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 5
            gradient1 = train_set'*(train_set*W(:,4)-train_label)/batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1)), gradient1, zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 1
            gradient1 = train_set'*(train_set*W(:,5)-train_label)/batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), zeros(size(gradient1)),zeros(size(gradient1)),gradient1];
        end
        gradientnew =  -opts.lambda2* invomega*W'*W*invomega   - muk * invomega;
        if i ==  1
            gradient = gradientnew;
        end
        gradient = 0.9*gradient + 0.1*gradientnew;

    elseif strcmp(opts.methods,'igt')
        Omegaz{1} = Omega{1} + (1 - gamma)/gamma*(Omega{1}-Omegaold{1});
        Wz = W + (1 - gamma)/gamma*(W-Wold);
        [U, S] = eig(Omegaz{1});
        invomega = U * diag(1./(diag(S)+0.01)) * U';
        tmpW = Wz*invomega;
        % invomega = inv(Omega);
        if mod(i,5) + 1 == 2
            gradient1 = train_set'*(train_set*W(:,1)-train_label)/batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1);
            gradient1 = [gradient1, zeros(size(gradient1)) ,zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 3
            gradient1 = train_set'*(train_set*W(:,2)-train_label)/batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2);
            gradient1 = [zeros(size(gradient1)), gradient1,  zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 4
            gradient1 = train_set'*(train_set*W(:,3)-train_label)/batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), gradient1,zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 5
            gradient1 = train_set'*(train_set*W(:,4)-train_label)/batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1)), gradient1, zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 1
            gradient1 = train_set'*(train_set*W(:,5)-train_label)/batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), zeros(size(gradient1)),zeros(size(gradient1)),gradient1];
        end
        gradient =  -opts.lambda2* invomega*Wz'*Wz*invomega   - muk * invomega;
    elseif strcmp(opts.methods,'recursiv')
        [U, S] = eig(Omega{1});
        invomega = U * diag(1./(diag(S)+0.01)) * U';
        tmpW = W*invomega;
        % invomega = inv(Omega);
        if mod(i,5) + 1 == 2
            gradient1 = train_set'*(train_set*W(:,1)-train_label)/batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1);
            gradient1 = [gradient1, zeros(size(gradient1)) ,zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 3
            gradient1 = train_set'*(train_set*W(:,2)-train_label)/batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2);
            gradient1 = [zeros(size(gradient1)), gradient1,  zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 4
            gradient1 = train_set'*(train_set*W(:,3)-train_label)/batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), gradient1,zeros(size(gradient1)),zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 5
            gradient1 = train_set'*(train_set*W(:,4)-train_label)/batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1)), gradient1, zeros(size(gradient1))];
        elseif mod(i,5) + 1 == 1
            gradient1 = train_set'*(train_set*W(:,5)-train_label)/batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5);
            gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), zeros(size(gradient1)),zeros(size(gradient1)),gradient1];
        end
        gradientnew =  -opts.lambda2* invomega*W'*W*invomega;
        if i ==  1
            gradientbar = gradientnew;
        end
        gradientbar = gradientnew + (1 - gamma)* (gradientbar - gradientnew);
        gradientbar = min(1, opts.radius)*gradientbar;
        gradient = gradientbar - muk * invomega;
    end

    tmp{1} = Omega{1}*gradient*Omega{1};
    lambda = -trans.Amap(tmp);
    schur = zeros(length(b),length(b));
    par.iter = i;

    schur = schurmat_sblk(blk,At,par,schur,1, Omega, Omega);


    if ~issparse(schur); schur = sparse(schur); end;
    L.matfct_options = 'spchol';
    [L.R,indef,L.perm] = chol(schur,'vector');
    L.Rt  = L.R';
    L.matdim = length(schur);
    diagR = full(diag(L.R)).^2;
    coeff.mat11 = schur;
    coeff.mat22 = [];
    coeff.mat12 = [];
    [lambda,resnrm,solve_ok] = symqmr(coeff,lambda,L,[],[],10);
    Omegaold = Omega;
    Wold = W;
    tmplambda = trans.ATmap(lambda);

    tmpgradient = Omega{1}*(gradient + tmplambda{1})*Omega{1};
    tmpall1 = data.train_set1'*(data.train_set1*W(:,1)-data.train_label1);
   tmpall1 = norm(2*tmpall1./(1+tmpall1.^2).^2/batchsize*1e4 + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1),'fro')^2 ;

        tmpall2 = data.train_set2'*(data.train_set2*W(:,2)-data.train_label2);
    tmpall2 = norm(2*tmpall2./(1+tmpall2.^2).^2/batchsize*1e4 + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1),'fro')^2 ;

        tmpall3 = data.train_set3'*(data.train_set3*W(:,3)-data.train_label3);
   tmpall3 = norm(2*tmpall3./(1+tmpall3.^2).^2/batchsize*1e4 + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1),'fro')^2 ;

        tmpall4 = data.train_set4'*(data.train_set4*W(:,4)-data.train_label3);
    tmpall4 = norm(2*tmpall4./(1+tmpall4.^2).^2/batchsize*1e4 + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1),'fro')^2 ;

        tmpall5 = data.train_set5'*(data.train_set5*W(:,5)-data.train_label4);
    tmpall5 = norm(2*tmpall5./(1+tmpall5.^2).^2/batchsize*1e4 + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1),'fro')^2 ;
    % tmpall1 = norm(data.train_set1'*(data.train_set1*W(:,1)-data.train_label1)/batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1),'fro')^2;
    % tmpall2 = norm(data.train_set2'*(data.train_set2*W(:,2)-data.train_label2)/batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2),'fro')^2;
    % tmpall3 = norm(data.train_set3'*(data.train_set3*W(:,3)-data.train_label3)/batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3),'fro')^2;
    % tmpall4 = norm(data.train_set4'*(data.train_set4*W(:,4)-data.train_label4)/batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4),'fro')^2;
    % tmpall5 = norm(data.train_set5'*(data.train_set5*W(:,5)-data.train_label5)/batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5),'fro')^2;

    tmperror = (norm(Omega{1}*(gradient + tmplambda{1}),'fro')^2 + tmpall1 + tmpall2 + tmpall3 + tmpall4 + tmpall5 )/(1 + norm(data.train_set1/batchsize/10,'fro') )

     f1 = data.train_set1*W(:,1) - data.train_label1;
    f2 = data.train_set2*W(:,2) - data.train_label2;
    f3 = data.train_set3*W(:,3) - data.train_label3;
    f4 = data.train_set4*W(:,4) - data.train_label4;
    f5 = data.train_set5*W(:,5) - data.train_label5;
    scale =1e4;
    tmpf1 = sum(f1.^2./(1 + f1.^2),'all')/batchsize*scale;
    tmpf2 = sum(f2.^2./(1 + f2.^2),'all')/batchsize*scale;
    tmpf3 = sum(f3.^2./(1 + f3.^2),'all')/batchsize*scale;
    tmpf4 = sum(f4.^2./(1 + f4.^2),'all')/batchsize*scale;
    tmpf5 = sum(f5.^2./(1 + f5.^2),'all')/batchsize*scale;
    % 
    fval = tmpf1 + tmpf2 + tmpf3 + tmpf4 + tmpf5 + opts.lambda1/2*norm(W,'fro')^2 + opts.lambda2/2*trace(W*Omega{1}*W')
    
    fval = norm(data.train_set1*W(:,1) - data.train_label1)^2/batchsize + norm(data.train_set2*W(:,2) - data.train_label2)^2/batchsize + norm(data.train_set3*W(:,3) - data.train_label3)^2/batchsize ...
        + norm(data.train_set4*W(:,4) - data.train_label4)^2/batchsize + + norm(data.train_set5*W(:,5) - data.train_label5)^2/batchsize + opts.lambda1/2*norm(W,'fro')^2 + opts.lambda2/2*trace(W*Omega{1}*W');
    % gradientave
    out.error = [out.error tmperror];
    out.fval = [out.fval fval];


    % tmpgradient = (tmpgradient + tmpgradient')/2;

    % x{1} = (x{1} + x{1}') /2;
    % if i <1000
    % lr = min(max(1/i,0.0001),1);
    % lr = min(max(6/i^(3/4),0.01),1);


    % else
    %     lr = max(1/i,0.001);
    % end

    if ~strcmp(opts.methods,'recursiv')
        % lr = min(max(6/i^(3/4),0.01),1);
        lr = opts.lr(i);
        lr2 = opts.lr(i);
        if tmperror<1e-3
            break;
        end
        Omega{1} = Omegaold{1} - opts.lr(i)*(tmpgradient)/(1e-6 + norm(Omega{1}*(gradient + tmplambda{1}),'fro'));
        Omega{1} = (Omega{1} + Omega{1}') /2;
        W = W-lr2*gradient1;
    else

        if tmperror<1e-4
            break;
        end
        lr = opts.lr(i)/(1e-6 +norm(Omega{1}*(gradient + tmplambda{1}),'fro'));
        lr2 = opts.lr(i);
        Omega{1} = Omegaold{1} - lr*tmpgradient;
        W = W-lr2*gradient1;
        Omega{1} = (Omega{1} + Omega{1}') /2;
    end


    muk = max(muk/1.2,1e-9);
end

out.Omega = Omega;
out.W = W;


end
% 
% function out = sto_ipm_multi_task(blk,x0,data,At,b,opts)
% par.numcolAt = length(b);
% % kappa = opts.kappa;
% % C{1} = Csto(:,:,1);
% W= x0;
% Omega{1} = opts.Omega0;
% z = Omega;
% 
% [At,~,Omega,z,par.permA,par.permZ] = sortA(blk,At,[],b,Omega,z);
% par.spdensity   = 0.4;
% par.smallblkdim = 50;
% [par.isspA,par.nzlistA,par.nzlistAsum,par.isspAy,par.nzlistAy] = nzlist(blk,At,par);
% batchsize = opts.batchsize;
% muk = 1;
% trans.Amap = opts.Amap;
% trans.ATmap = opts.ATmap;
% gamma = opts.gamma;
% lr = opts.lr;
% x = x0;
% Omegaold = Omega;
% Wold = W;
% out.error = [];
% out.fval = [];
% % gradientave = sum(Csto(:,:,:),3)/batchsize;
% for i = 1:opts.maxiter
%     if mod(i,5) + 1 == 2
%         train_set = data.train_set1;
%         train_label = data.train_label1;
%         batchsize = length(train_label);
%     elseif mod(i,5) + 1 == 3
%         train_set = data.train_set2;
%         train_label = data.train_label2;
%         batchsize = length(train_label);
%     elseif mod(i,5) + 1 == 4
%         train_set = data.train_set3;
%         train_label = data.train_label3;
%         batchsize = length(train_label);
%     elseif mod(i,5) + 1 == 5
%         train_set = data.train_set4;
%         train_label = data.train_label4;
%         batchsize = length(train_label);
%     elseif mod(i,5) + 1 == 1
%         train_set = data.train_set5;
%         train_label = data.train_label5;
%         batchsize = length(train_label);
%     end
%     % batchnum = mod(i,10);
%     % batchnum = 0;
%     % batch = batchnum*batchsize + 1: batchsize*(batchnum+1);
%     if strcmp(opts.methods,'sto')
%         [U, S] = eig(Omega{1});
%         invomega = U * diag(1./(diag(S)+0.01)) * U';
%         tmpW = W*invomega;
%         % invomega = inv(Omega);
%         if mod(i,5) + 1 == 2
%             tmpg = train_set'*(train_set*W(:,1)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2  + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1);
%             gradient1 = [gradient1, zeros(size(gradient1)) ,zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 3
%             tmpg = train_set'*(train_set*W(:,2)-train_label);
%             gradient1 =2*tmpg.^2./(1+tmpg.^2).^2  + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2);
%             gradient1 = [zeros(size(gradient1)), gradient1,  zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 4
%             tmpg = train_set'*(train_set*W(:,3)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2  + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), gradient1,zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 5
%             tmpg = train_set'*(train_set*W(:,4)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2+ opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1)), gradient1, zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 1
%             tmpg = train_set'*(train_set*W(:,5)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2  + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), zeros(size(gradient1)),zeros(size(gradient1)),gradient1];
%         end
%         gradient =  -opts.lambda2* invomega*W'*W*invomega   - muk * invomega;
%     elseif strcmp(opts.methods,'mom')
%         [U, S] = eig(Omega{1});
%         invomega = U * diag(1./(diag(S)+0.01)) * U';
%         tmpW = W*invomega;
%         % invomega = inv(Omega);
%         if mod(i,5) + 1 == 2
%             tmpg = train_set'*(train_set*W(:,1)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1);
%             gradient1 = [gradient1, zeros(size(gradient1)) ,zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 3
%             tmpg = train_set'*(train_set*W(:,2)-train_label);
%             gradient1 =2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2);
%             gradient1 = [zeros(size(gradient1)), gradient1,  zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 4
%             tmpg = train_set'*(train_set*W(:,3)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), gradient1,zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 5
%             tmpg = train_set'*(train_set*W(:,4)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1)), gradient1, zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 1
%             tmpg = train_set'*(train_set*W(:,5)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), zeros(size(gradient1)),zeros(size(gradient1)),gradient1];
%         end
%         gradientnew =  -opts.lambda2* invomega*W'*W*invomega   - muk * invomega;
%         if i ==  1
%             gradient = gradientnew;
%         end
%         gradient = 0.9*gradient + 0.1*gradientnew;
% 
%     elseif strcmp(opts.methods,'igt')
%         Omegaz{1} = Omega{1} + (1 - gamma)/gamma*(Omega{1}-Omegaold{1});
%         Wz = W + (1 - gamma)/gamma*(W-Wold);
%         [U, S] = eig(Omegaz{1});
%         invomega = U * diag(1./(diag(S)+0.01)) * U';
%         tmpW = Wz*invomega;
%         % invomega = inv(Omega);
%         if mod(i,5) + 1 == 2
%             tmpg = train_set'*(train_set*W(:,1)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1);
%             gradient1 = [gradient1, zeros(size(gradient1)) ,zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 3
%             tmpg = train_set'*(train_set*W(:,2)-train_label);
%             gradient1 =2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2);
%             gradient1 = [zeros(size(gradient1)), gradient1,  zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 4
%             tmpg = train_set'*(train_set*W(:,3)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), gradient1,zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 5
%             tmpg = train_set'*(train_set*W(:,4)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1)), gradient1, zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 1
%             tmpg = train_set'*(train_set*W(:,5)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), zeros(size(gradient1)),zeros(size(gradient1)),gradient1];
%         end
%         gradient =  -opts.lambda2* invomega*Wz'*Wz*invomega   - muk * invomega;
%     elseif strcmp(opts.methods,'recursiv')
%         [U, S] = eig(Omega{1});
%         invomega = U * diag(1./(diag(S)+0.01)) * U';
%         tmpW = W*invomega;
%         % invomega = inv(Omega);
%         if mod(i,5) + 1 == 2
%             tmpg = train_set'*(train_set*W(:,1)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1);
%             gradient1 = [gradient1, zeros(size(gradient1)) ,zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 3
%             tmpg = train_set'*(train_set*W(:,2)-train_label);
%             gradient1 =2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2);
%             gradient1 = [zeros(size(gradient1)), gradient1,  zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 4
%             tmpg = train_set'*(train_set*W(:,3)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), gradient1,zeros(size(gradient1)),zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 5
%             tmpg = train_set'*(train_set*W(:,4)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)),zeros(size(gradient1)), gradient1, zeros(size(gradient1))];
%         elseif mod(i,5) + 1 == 1
%             tmpg = train_set'*(train_set*W(:,5)-train_label);
%             gradient1 = 2*tmpg.^2./(1+tmpg.^2).^2 /batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5);
%             gradient1 = [zeros(size(gradient1)),zeros(size(gradient1)), zeros(size(gradient1)),zeros(size(gradient1)),gradient1];
%         end
%         gradientnew =  -opts.lambda2* invomega*W'*W*invomega;
%         if i ==  1
%             gradientbar = gradientnew;
%         end
%         gradientbar = gradientnew + (1 - gamma)* (gradientbar - gradientnew);
%         gradientbar = min(1, opts.radius)*gradientbar;
%         gradient = gradientbar - muk * invomega;
%     end
% 
%     tmp{1} = Omega{1}*gradient*Omega{1};
%     lambda = -trans.Amap(tmp);
%     schur = zeros(length(b),length(b));
%     par.iter = i;
% 
%     schur = schurmat_sblk(blk,At,par,schur,1, Omega, Omega);
% 
% 
%     if ~issparse(schur); schur = sparse(schur); end;
%     L.matfct_options = 'spchol';
%     [L.R,indef,L.perm] = chol(schur,'vector');
%     L.Rt  = L.R';
%     L.matdim = length(schur);
%     diagR = full(diag(L.R)).^2;
%     coeff.mat11 = schur;
%     coeff.mat22 = [];
%     coeff.mat12 = [];
%     [lambda,resnrm,solve_ok] = symqmr(coeff,lambda,L,[],[],10);
%     Omegaold = Omega;
%     Wold = W;
%     tmplambda = trans.ATmap(lambda);
% 
%     tmpgradient = Omega{1}*(gradient + tmplambda{1})*Omega{1};
%     tmpf1 = data.train_set1'*(data.train_set1*W(:,1)-data.train_label1);
%     tmpf2 = data.train_set2'*(data.train_set1*W(:,2)-data.train_label2);
%     tmpf3 = data.train_set3'*(data.train_set1*W(:,3)-data.train_label3);
%     tmpf4 = data.train_set4'*(data.train_set1*W(:,4)-data.train_label4);
%     tmpf5 = data.train_set5'*(data.train_set1*W(:,5)-data.train_label5);
% 
%     tmpg1  =  2*tmpf1./(1+tmpf1.^2).^2;
%     tmpg2  =  2*tmpf2./(1+tmpf2.^2).^2;
%     tmpg3  =  2*tmpf3./(1+tmpf3.^2).^2;
%     tmpg4  =  2*tmpf4./(1+tmpf4.^2).^2;
%     tmpg5  =  2*tmpf5./(1+tmpf5.^2).^2;
% 
%     tmpall1 = norm(tmpf1/batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1),'fro')^2;
%     tmpall2 = norm(tmpf2/batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2),'fro')^2;
%     tmpall3 = norm(tmpf3/batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3),'fro')^2;
%     tmpall4 = norm(tmpf4/batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4),'fro')^2;
%     tmpall5 = norm(tmpf5/batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5),'fro')^2;
% 
%     tmperror = (norm(Omega{1}*(gradient + tmplambda{1}),'fro')^2 + tmpall1 + tmpall2 + tmpall3 + tmpall4 + tmpall5 )/(1 + norm(data.train_set1/batchsize/10,'fro') )
% 
%     f1 = data.train_set1*W(:,1) - data.train_label1;
%     f2 = data.train_set2*W(:,2) - data.train_label2;
%     f3 = data.train_set3*W(:,3) - data.train_label3;
%     f4 = data.train_set4*W(:,4) - data.train_label4;
%     f5 = data.train_set5*W(:,5) - data.train_label5;
% 
%     f1 = sum(f1.^2/(1 + f1.^2),'all')/batchsize;
%     f2 = sum(f2.^2/(1 + f2.^2),'all')/batchsize;
%     f3 = sum(f3.^2/(1 + f3.^2),'all')/batchsize;
%     f4 = sum(f4.^2/(1 + f4.^2),'all')/batchsize;
%     f5 = sum(f5.^2/(1 + f5.^2),'all')/batchsize;
% 
%     fval = f1 + f2 + f3 + f4 + f5 + opts.lambda1/2*norm(W,'fro')^2 + opts.lambda2/2*trace(W*Omega{1}*W');
%     % gradientave
%     out.error = [out.error tmperror];
%     out.fval = [out.fval fval];
% 
% 
%     % tmpgradient = (tmpgradient + tmpgradient')/2;
% 
%     % x{1} = (x{1} + x{1}') /2;
%     % if i <1000
%     % lr = min(max(1/i,0.0001),1);
%     % lr = min(max(6/i^(3/4),0.01),1);
% 
% 
%     % else
%     %     lr = max(1/i,0.001);
%     % end
% 
%     if ~strcmp(opts.methods,'recursiv')
%         % lr = min(max(6/i^(3/4),0.01),1);
%         lr = opts.lr(i);
%         lr2 = opts.lr(i);
%         if tmperror<1e-6
%             break;
%         end
%         Omega{1} = Omegaold{1} - opts.lr(i)*(tmpgradient)/(1e-6 + norm(Omega{1}*(gradient + tmplambda{1}),'fro'));
%         Omega{1} = (Omega{1} + Omega{1}') /2;
%         W = W-lr2*gradient1;
%     else
% 
%         if tmperror<1e-6
%             break;
%         end
%         lr = opts.lr(i)/(1e-6 +norm(Omega{1}*(gradient + tmplambda{1}),'fro'));
%         lr2 = opts.lr(i);
%         Omega{1} = Omegaold{1} - lr*tmpgradient;
%         W = W-lr2*gradient1;
%         Omega{1} = (Omega{1} + Omega{1}') /2;
%     end
% 
% 
%     muk = max(muk/1.2,1e-9);
% end
% 
% out.x = x;
% 
% 
% end
