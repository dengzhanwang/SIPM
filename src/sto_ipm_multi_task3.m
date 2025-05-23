function out = sto_ipm_multi_task3(blk,x0,data,At,b,opts)
par.numcolAt = length(b);

W= x0;
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
trainsetall = data.trainsetall;
trainlabelall = data.trainlabelall;
numall = opts.numall;

for i = 1:opts.maxiter
     if opts.options == 1
    batchnum = mod(i,numall)+1;
    % batchnum = 0;
    batch = (batchnum-1)*batchsize + 1: batchsize*(batchnum);
     elseif opts.options == 2
        batch = 1:min(i*10,size(trainsetall,1));
        batchnum = length(batch);
     elseif opts.options == 3
        batch = 1:size(trainsetall,1);
        batchnum = length(batch);
     end

    if strcmp(opts.methods,'sto')
        [U, S] = eig(Omega{1});
        invomega = U * diag(1./(diag(S)+0.01)) * U';
        tmpW = W*invomega;
        if opts.options == 1
        gradient1tmp = trainsetall(batch,:)'*(trainsetall(batch,:)*W(:,batchnum)-trainlabelall(batch,:)); 
        gradient1tmp = 2*gradient1tmp./(1+gradient1tmp.^2).^2/batchsize + opts.lambda1*W(:,batchnum) + opts.lambda2*tmpW(:,batchnum);
        gradient1 = zeros(size(W)); 
        gradient1(:,batchnum) = gradient1tmp;
        else
            gradient1 = zeros(size(W)); 
        for jj = 1:ceil(batchnum/batchsize)
            
            if jj == ceil(i/batchsize)
                    batch = (jj-1)*batchsize + 1: batchnum;
                    gradient1tmp = trainsetall(batch,:)'*(trainsetall(batch,:)*W(:,jj)-trainlabelall(batch,:)); 
            else
                batch = (jj-1)*batchsize + 1: batchsize*(jj);
            gradient1tmp = trainsetall(batch,:)'*(trainsetall(batch,:)*W(:,jj)-trainlabelall(batch,:)); 
            end
        gradient1tmp = 2*gradient1tmp./(1+gradient1tmp.^2).^2/batchsize + opts.lambda1*W(:,jj) + opts.lambda2*tmpW(:,jj);
        
        gradient1(:,jj) = gradient1tmp;
        end
        end

        gradient =  -opts.lambda2* invomega*W'*W*invomega   - muk * invomega;
    elseif strcmp(opts.methods,'mom')
        [U, S] = eig(Omega{1});
        invomega = U * diag(1./(diag(S)+0.01)) * U';
        tmpW = W*invomega;
        gradient1tmp = trainsetall(batch,:)'*(trainsetall(batch,:)*W(:,batchnum)-trainlabelall(batch,:)); 
        gradient1tmp = 2*gradient1tmp./(1+gradient1tmp.^2).^2/batchsize + opts.lambda1*W(:,batchnum) + opts.lambda2*tmpW(:,batchnum);
        gradient1 = zeros(size(W)); 
        gradient1(:,batchnum) = gradient1tmp;
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
           gradient1tmp = trainsetall(batch,:)'*(trainsetall(batch,:)*W(:,batchnum)-trainlabelall(batch,:)); 
        gradient1tmp = 2*gradient1tmp./(1+gradient1tmp.^2).^2/batchsize + opts.lambda1*W(:,batchnum) + opts.lambda2*tmpW(:,batchnum);
        gradient1 = zeros(size(W)); 
        gradient1(:,batchnum) = gradient1tmp;
        gradient =  -opts.lambda2* invomega*Wz'*Wz*invomega   - muk * invomega;
    elseif strcmp(opts.methods,'recursiv')
        [U, S] = eig(Omega{1});
        invomega = U * diag(1./(diag(S)+0.01)) * U';
        tmpW = W*invomega;
            gradient1tmp = trainsetall(batch,:)'*(trainsetall(batch,:)*W(:,batchnum)-trainlabelall(batch,:)); 
        gradient1tmp = 2*gradient1tmp./(1+gradient1tmp.^2).^2/batchsize + opts.lambda1*W(:,batchnum) + opts.lambda2*tmpW(:,batchnum);
        gradient1 = zeros(size(W)); 
        gradient1(:,batchnum) = gradient1tmp;
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
    tmpg1 = data.train_set1'*(data.train_set1*W(:,1)-data.train_label1);
    tmpg1 = 2*tmpg1./(1 + tmpg1.^2).^2;
    tmpg2 = data.train_set1'*(data.train_set2*W(:,2)-data.train_label2);
    tmpg2 = 2*tmpg2./(1 + tmpg2.^2).^2;
    tmpg3 = data.train_set1'*(data.train_set3*W(:,3)-data.train_label3);
    tmpg3 = 2*tmpg3./(1 + tmpg3.^2).^2;
    tmpg4 = data.train_set1'*(data.train_set4*W(:,4)-data.train_label4);
    tmpg4 = 2*tmpg4./(1 + tmpg4.^2).^2;
    tmpg5 = data.train_set1'*(data.train_set5*W(:,5)-data.train_label5);
    tmpg5 = 2*tmpg5./(1 + tmpg5.^2).^2;
    tmpall1 = norm(tmpg1/batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1),'fro')^2;
    tmpall2 = norm(tmpg2/batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2),'fro')^2;
    tmpall3 = norm(tmpg3/batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3),'fro')^2;
    tmpall4 = norm(tmpg4/batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4),'fro')^2;
    tmpall5 = norm(tmpg5/batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5),'fro')^2;
    % if i == 1
    %     tmperror0 = (norm(Omega{1}*(gradient + tmplambda{1}),'fro')^2 + tmpall1 + tmpall2 + tmpall3 + tmpall4 + tmpall5 );
    % end
    tmperror = sqrt(norm(Omega{1}*(gradient + tmplambda{1}),'fro')^2 + tmpall1 + tmpall2 + tmpall3 + tmpall4 + tmpall5 );

     f1 = data.train_set1*W(:,1) - data.train_label1;
    f2 = data.train_set2*W(:,2) - data.train_label2;
    f3 = data.train_set3*W(:,3) - data.train_label3;
    f4 = data.train_set4*W(:,4) - data.train_label4;
    f5 = data.train_set5*W(:,5) - data.train_label5;

    tmpf1 = sum(f1.^2./(1 + f1.^2),'all')/batchsize;
    tmpf2 = sum(f2.^2./(1 + f2.^2),'all')/batchsize;
    tmpf3 = sum(f3.^2./(1 + f3.^2),'all')/batchsize;
    tmpf4 = sum(f4.^2./(1 + f4.^2),'all')/batchsize;
    tmpf5 = sum(f5.^2./(1 + f5.^2),'all')/batchsize;
    % 
    fval = tmpf1*2 + tmpf2*2 + tmpf3*2 + tmpf4*2 + tmpf5*2 + opts.lambda1/2*norm(W,'fro')^2 + opts.lambda2/2*trace(W*Omega{1}*W');
    
    % fval = norm(data.train_set1*W(:,1) - data.train_label1)^2/batchsize + norm(data.train_set2*W(:,2) - data.train_label2)^2/batchsize + norm(data.train_set3*W(:,3) - data.train_label3)^2/batchsize ...
        % + norm(data.train_set4*W(:,4) - data.train_label4)^2/batchsize + + norm(data.train_set5*W(:,5) - data.train_label5)^2/batchsize + opts.lambda1/2*norm(W,'fro')^2 + opts.lambda2/2*trace(W*Omega{1}*W');
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
        % if tmperror<1e-3
        %     break;
        % end
        Omega{1} = Omegaold{1} - opts.lr(i)*(tmpgradient)/(1e-6 + norm(Omega{1}*(gradient + tmplambda{1}),'fro'));
        Omega{1} = (Omega{1} + Omega{1}') /2;
        W = W-lr2*gradient1*2;
    else

        % if tmperror<1e-4
        %     break;
        % end
        lr = opts.lr(i)/(1e-6 +norm(Omega{1}*(gradient + tmplambda{1}),'fro'));
        lr2 = opts.lr(i);
        Omega{1} = Omegaold{1} - lr*tmpgradient;
        W = W-lr2*gradient1*2;
        Omega{1} = (Omega{1} + Omega{1}') /2;
    end


    muk = max(muk/1.2,1e-9);
end

out.Omega = Omega;
out.W = W;


end
% 