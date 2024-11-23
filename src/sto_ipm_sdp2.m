function out = sto_ipm_sdp2(blk,x0,At,Csto,b,opts)

par.numcolAt = length(b);
kappa = opts.kappa;
C{1} = Csto(:,:,1);
x= x0;
z = x0;
[At,C,x,z,par.permA,par.permZ] = sortA(blk,At,C,b,x,z);
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
xold = x;
for i = 1:opts.maxiter
    % if strcmp(blk{1},'s')
    %     invxv = inv(x);
    % elseif strcmp(blk{1},'l')
    %     invxv = 1./x;
    % end
    batch = opts.batch(1): opts.batch(1) + batchsize-1;
    if strcmp(opts.methods,'sto')
        [U, S] = eig(x{1});
        invxv = U * diag(1./diag(S)) * U';
        gradient = sum(Csto(:,:,batch),3)/batchsize + kappa* U * diag(diag(S)./diag(1+S)) * U'   + muk * invxv;
    elseif strcmp(opts.methods,'mom')
        [U, S] = eig(x{1});
        invxv = U * diag(1./diag(S)) * U';
        gradientnew = sum(Csto(:,:,batch),3)/batchsize + kappa* U * diag(diag(S)./diag(1+S)) * U'   + muk * invxv;
        if i ==  1
            gradient = gradientnew;
        end
        gradient = 0.9*gradient + 0.1*gradientnew;
    elseif strcmp(opts.methods,'igt')
        z{1} = x{1} + (1 - gamma)/gamma*(x{1}-xold{1});
        [U, S] = eig(z{1});
        invxv = U * diag(1./diag(S)) * U';
        gradientnew = sum(Csto(:,:,batch),3)/batchsize + kappa* U * diag(diag(S)./diag(1+S)) * U'   + muk * invxv;
        if i ==  1
            gradient = gradientnew;
        end
        gradient = (1-gamma)*gradient + gamma*gradientnew;
    elseif strcmp(opts.methods,'recursiv')
        [U, S] = eig(x{1});
        invxv = U * diag(1./diag(S)) * U';
        gradientnew = sum(Csto(:,:,batch),3)/batchsize + kappa* U * diag(diag(S)./diag(1+S)) * U' ;
        if i ==  1
            gradientbar = gradientnew;
        end
        gradientbar = gradientnew + (1 - gamma)* (gradientbar - gradientnew);
        gradientbar = min(1, opts.radius/norm(gradientbar*x{1},'fro'))*gradientbar;
        gradient = gradientbar - muk * invxv;
    end

    if strcmp(blk{1},'s')
        tmp{1} = x{1}*gradient*x{1};
        lambda = -trans.Amap(tmp);
    elseif strcmp(blk{1},'l')
        lambda = -trans.Amap(gradient.*x.^2);
    end

    schur = zeros(length(b),length(b));
    par.iter = i;
    schur = schurmat_sblk(blk,At,par,schur,1,x,x);
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
    xold = x;
    tmplambda = trans.ATmap(lambda);
    tmpgradient = x{1}*(gradient + tmplambda{1})*x{1};
    % tmpgradient = (tmpgradient + tmpgradient')/2;
    norm(x{1}*(gradient + tmplambda{1}),'fro')/(1 + norm(sum(Csto(:,:,batch),3)/batchsize,'fro') )
    % x{1} = (x{1} + x{1}') /2;
    if ~strcmp(opts.methods,'recursiv')
        x{1} = xold{1} - lr*(tmpgradient)/norm(x{1}*(gradient + tmplambda{1}),'fro');
    else
        x{1} = xold{1} - lr*tmpgradient;
    end
    x{1} = (x{1} + x{1}') /2;
    if i <500
        lr = min(max(6/sqrt(i),0.01),1);
    else
        lr = max(1/i,0.001);
    end
    muk = max(muk/1.1,1e-9);
end


end