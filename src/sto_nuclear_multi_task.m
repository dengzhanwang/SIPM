function out = sto_ipm_multi_task(blk,x0,data,At,b,opts)
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
numall = size(data.train_set1,1)/batchsize;
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
    if strcmp(opts.methods,'sto')
        [U, S] = eig(Omega{1});
        invomega = U * diag(1./(diag(S)+0.01)) * U';
        tmpW = W*invomega;
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
    tmpall1 = norm(data.train_set1'*(data.train_set1*W(:,1)-data.train_label1)/batchsize + opts.lambda1*W(:,1) + opts.lambda2*tmpW(:,1),'fro')^2;
    tmpall2 = norm(data.train_set2'*(data.train_set2*W(:,2)-data.train_label2)/batchsize + opts.lambda1*W(:,2) + opts.lambda2*tmpW(:,2),'fro')^2;
    tmpall3 = norm(data.train_set3'*(data.train_set3*W(:,3)-data.train_label3)/batchsize + opts.lambda1*W(:,3) + opts.lambda2*tmpW(:,3),'fro')^2;
    tmpall4 = norm(data.train_set4'*(data.train_set4*W(:,4)-data.train_label4)/batchsize + opts.lambda1*W(:,4) + opts.lambda2*tmpW(:,4),'fro')^2;
    tmpall5 = norm(data.train_set5'*(data.train_set5*W(:,5)-data.train_label5)/batchsize + opts.lambda1*W(:,5) + opts.lambda2*tmpW(:,5),'fro')^2;

    tmperror = (norm(Omega{1}*(gradient + tmplambda{1}),'fro')^2 + tmpall1 + tmpall2 + tmpall3 + tmpall4 + tmpall5 )/(1 + norm(data.train_set1/batchsize/10,'fro') );

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
    fval = tmpf1 + tmpf2 + tmpf3 + tmpf4 + tmpf5 + opts.lambda1/2*norm(W,'fro')^2 + opts.lambda2/2*trace(W*Omega{1}*W');
    
    out.error = [out.error tmperror];
    out.fval = [out.fval fval];






    if ~strcmp(opts.methods,'recursiv')
        lr = opts.lr(i);
        lr2 = opts.lr(i);

        W = W-lr2*gradient1;
    else
        lr = opts.lr(i)/(1e-6 +norm(Omega{1}*(gradient + tmplambda{1}),'fro'));
        lr2 = opts.lr(i);
        W = W-lr2*gradient1;
    end


    muk = max(muk/1.2,1e-9);
end

out.Omega = Omega;
out.W = W;


end
