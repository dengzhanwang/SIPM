function out = sto_ipmsgd_full(blk,x0,At,Csto,b,opts)
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
out.error = [];
out.fval = [];
out.outiter = opts.maxiter;
tmperrorsum = 0
% gradientave = sum(Csto(:,:,:),3)/batchsize;
for i = 1:opts.maxiter
    batchnum = mod(i,10);
    % batchnum = 0;
    batch = batchnum*batchsize + 1: batchsize*(batchnum+1);
    if strcmp(opts.methods,'sto')
        if strcmp(blk{1},'s')
            [U, S] = eig(x{1});
            invxv = U * diag(1./diag(S)) * U';
            gradient = sum(Csto(:,:,batch),3)/batchsize + kappa* U * diag(diag(S)./diag(1+S)) * U'  - muk * invxv;
        elseif strcmp(blk{1},'l')
            tmpg = [ ones(length(x{1})/2,1) ; zeros(length(x{1})/2,1) ];
            gradient = sum(Csto(:,:,batch),3)/batchsize*x{1} - tmpg  - muk * 1./x{1};
        end
    elseif strcmp(opts.methods,'mom')
        if strcmp(blk{1},'s')
            [U, S] = eig(x{1});
            invxv = U * diag(1./diag(S)) * U';
            gradientnew = sum(Csto(:,:,batch),3)/batchsize + kappa* U * diag(diag(S)./diag(1+S)) * U'   - muk * invxv;
        elseif strcmp(blk{1},'l')
            tmpg = [ ones(length(x{1})/2,1) ; zeros(length(x{1})/2,1) ];
            gradientnew = sum(Csto(:,:,batch),3)/batchsize*x{1} - tmpg  - muk * 1./x{1};
        end
        if i ==  1
            gradient = gradientnew;
        end
        gradient = 0.9*gradient + 0.1*gradientnew;
    elseif strcmp(opts.methods,'igt')
        z{1} = x{1} + (1 - gamma)/gamma*(x{1}-xold{1});
        if strcmp(blk{1},'s')
            [U, S] = eig(z{1});
            invxv = U * diag(1./diag(S)) * U';
            gradientnew = sum(Csto(:,:,batch),3)/batchsize + kappa* U * diag(diag(S)./diag(1+S)) * U'   - muk * invxv;
        elseif strcmp(blk{1},'l')
            tmpg = [ ones(length(x{1})/2,1) ; zeros(length(x{1})/2,1) ];
            gradientnew = sum(Csto(:,:,batch),3)/batchsize*z{1} - tmpg  - muk * 1./z{1};
        end
        if i ==  1
            gradient = gradientnew;
        end
        gradient = (1-gamma)*gradient + gamma*gradientnew;
    elseif strcmp(opts.methods,'recursiv')
        if strcmp(blk{1},'s')
            [U, S] = eig(x{1});
            invxv = U * diag(1./diag(S)) * U';
            gradientnew = sum(Csto(:,:,batch),3)/batchsize +  kappa* U * diag(diag(S)./diag(1+S)) * U'  ;
            if i ==  1
                gradientbar = gradientnew;
            end
            gradientbar = gradientnew + (1 - gamma)* (gradientbar - gradientnew);
            gradientbar = min(1, opts.radius/norm(gradientbar*x{1},'fro'))*gradientbar;
            gradient = gradientbar - muk * invxv;
        elseif strcmp(blk{1},'l')
            tmpg = [ ones(length(x{1})/2,1) ; zeros(length(x{1})/2,1) ];
            gradientnew = sum(Csto(:,:,batch),3)/batchsize*x{1} - tmpg;
            if i ==  1
                gradientbar = gradientnew;
            end
            gradientbar = gradientnew + (1 - gamma)* (gradientbar - gradientnew);
            gradientbar = min(1, opts.radius/norm(gradientbar.*x{1},'fro'))*gradientbar;
            gradient = gradientbar -  muk * 1./x{1};
        end

    end

    if strcmp(blk{1},'s')
        tmp{1} = x{1}*gradient*x{1};
        lambda = -trans.Amap(tmp);
    elseif strcmp(blk{1},'l')
        tmp{1} =  gradient.*x{1}.^2;
        lambda = -trans.Amap(tmp);
    elseif strcmp(blk{1},'q')
    end

    schur = zeros(length(b),length(b));
    par.iter = i;
    if strcmp(blk{1},'s')
        schur = schurmat_sblk(blk,At,par,schur,1,x,x);
    elseif strcmp(blk{1},'l')
        xtmp{1} = x{1}.^2;
        schur = schurmat_lblk(blk,At,par,schur,[],[],1,xtmp);
    end
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
    if strcmp(blk{1},'s')
        tmpgradient = x{1}*(gradient + tmplambda{1})*x{1};
        tmperrornew = norm(x{1}*(gradient + tmplambda{1}),'fro')^2;
        if batchnum~= 0 
        tmperrorsum = tmperrorsum + tmperrornew;
        tmperror = sqrt(tmperrorsum/batchnum)/(1 + norm(sum(Csto(:,:,:),3)/batchsize/10,'fro') ) 
        else
        tmperrorsum = tmperrorsum + tmperrornew;
        tmperror = sqrt(tmperrorsum/10)/(1 + norm(sum(Csto(:,:,:),3)/batchsize/10,'fro') ) 
        tmperrorsum  = 0;
        end
        fval = sum(sum(Csto(:,:,:),3).*x{1},'all');
        % gradientave
        out.error = [out.error tmperror];
        out.fval = [out.fval fval];
    elseif strcmp(blk{1},'l')
        tmpgradient = (gradient + tmplambda{1}).*x{1}.^2;
        tmperror = norm(x{1}.*(gradient + tmplambda{1}),'fro')/(1 + norm(sum(Csto(:,:,:),3)/batchsize,'fro') )
        fval = -sum(x{1}(1:length(x{1})/2),"all" ) + 0.5*x{1}'*sum(Csto(:,:,:),3)*x{1}/batchsize
        out.error = [out.error tmperror];
        out.fval = [out.fval fval];
    end


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
        % lr = opts.lr(i);
        if tmperror<1e-5
            break;
        end
        if strcmp(blk{1},'s')
            x{1} = xold{1} - opts.lr(i)*(tmpgradient)/norm(x{1}*(gradient + tmplambda{1}),'fro');
            x{1} = (x{1} + x{1}') /2;
        elseif strcmp(blk{1},'l')
            x{1} = xold{1} - opts.lr(i)*(tmpgradient)/norm(x{1}.*(gradient + tmplambda{1}),'fro');
        end
    else
        
        if tmperror<1e-5
            out.outiter = i;
        end

       if tmperror<1e-6
            break;
        end
        if strcmp(blk{1},'s')
            lr = opts.lr(i)/norm(x{1}*(gradient + tmplambda{1}),'fro');
            x{1} = xold{1} - lr*tmpgradient;
            x{1} = (x{1} + x{1}') /2;
        elseif strcmp(blk{1},'l')
            lr = opts.lr(i)/norm(x{1}.*(gradient + tmplambda{1}),'fro');
            x{1} = xold{1} - lr*tmpgradient;
        end
    end


    muk = max(muk/1.2,1e-9);
end

out.x = x;


end