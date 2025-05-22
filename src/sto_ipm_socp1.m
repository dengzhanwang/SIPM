function out = sto_ipm_socp1(blk,x0,At,data,opts)
lambda1 = opts.lambda1;
lambda2 = opts.lambda2;
z = x0;
n = opts.n;
batchsize = opts.batchsize;
muk = 1;

gamma = opts.gamma;
% % lr = opts.lr;
x = x0;
xold = x;
out.error = [];
out.fval = [];
tmperrorsum = 0;
numall = size(data.train_miss,1)/batchsize;
for i = 1:opts.maxiter
    if opts.options == 1
    batchnum = mod(i,numall);
    batch = batchnum*batchsize + 1: batchsize*(batchnum+1);
    elseif opts.options == 2
        batch = 1:min(i,size(data.train_miss,1));
        batchnum = 0;
    elseif opts.options == 3
        batch = 1:size(data.train_miss,1);
        batchnum = 0;        
    end
    if strcmp(opts.methods,'sto')
        det1 = x(1)^2 - x(2:n+1)'*x(2:n+1);
        tmp1 = [x(1); -x(2:n+1)]/det1;
        det2 = x(n+2)^2 - x(n+3:2*n+2)'*x(n+3:2*n+2);

        tmp2 = [x(n+2); -x(n+3:2*n+2)]/det2;
        invx = [tmp1; tmp2];
        tmpgradient = data.train_miss(batch,:)'*(data.train_miss(batch,:)*x(2:n+1)-data.train_label(batch))/length(batch);
        tmpgradient = 2*tmpgradient./(1+tmpgradient.^2).^2;
        tmpgradient = [0; tmpgradient; zeros(size(tmpgradient,1)+1,1) ];
        lambdagradient1 = zeros(2*n+2,1);
        lambdagradient1(1) = lambda1;
        lambdagradient1(n+2) = lambda2;
        gradient = tmpgradient +  lambdagradient1 - muk * (invx + tmpgradient);
    elseif strcmp(opts.methods,'mom')
        % x(1) = min( norm(x(2:n+1))+1,x(1) );
        det1 = x(1)^2 - x(2:n+1)'*x(2:n+1);
        tmp1 = [x(1); -x(2:n+1)]/det1;
        % x(n+2) = min( norm(x(n+3:2*n+2))+1,x(n+2) );
        det2 = x(n+2)^2 - x(n+3:2*n+2)'*x(n+3:2*n+2);
        tmp2 = [x(n+2); -x(n+3:2*n+2)]/det2;
        invx = [tmp1; tmp2];
        tmpgradient = data.train_miss(batch,:)'*(data.train_miss(batch,:)*x(2:n+1)-data.train_label(batch))/batchsize;
        tmpgradient = 2*tmpgradient./(1+tmpgradient.^2).^2;
        tmpgradient = [0; tmpgradient; zeros(size(tmpgradient,1)+1,1) ];
        % tmpgradient = tmpC'*(tmpC -y) /batchsize;
        lambdagradient1 = zeros(2*n+2,1);
        % lambdagradient2 = zeros(2*n,1);
        lambdagradient1(1) = lambda1;
        lambdagradient1(n+2) = lambda2;
        gradientnew = tmpgradient +  lambdagradient1 - muk * (invx + tmpgradient);
        if i ==  1
            gradient = gradientnew;
        end
        gradient = 0.9*gradient + 0.1*gradientnew;
    elseif strcmp(opts.methods,'igt')
        z = x + (1 - gamma)/gamma*(x-xold);
        det1 = z(1)^2 - z(2:n+1)'*z(2:n+1);
        tmp1 = [z(1); -z(2:n+1)]/det1;
        % x(n+2) = min( norm(x(n+3:2*n+2))+1,x(n+2) );
        det2 = z(n+2)^2 - z(n+3:2*n+2)'*z(n+3:2*n+2);
        tmp2 = [z(n+2); -z(n+3:2*n+2)]/det2;
        invx = [tmp1; tmp2];
        tmpgradient = data.train_miss(batch,:)'*(data.train_miss(batch,:)*x(2:n+1)-data.train_label(batch))/batchsize;
        tmpgradient = 2*tmpgradient./(1+tmpgradient.^2).^2;
        tmpgradient = [0; tmpgradient; zeros(size(tmpgradient,1)+1,1) ];
        tt = tmpgradient;
        lambdagradient1 = zeros(2*n+2,1);
        lambdagradient1(1) = lambda1;
        lambdagradient1(n+2) = lambda2;
        gradientnew = tmpgradient +  lambdagradient1 - muk * (invx + tmpgradient);
        if i ==  1
            gradient = gradientnew;
        end
        gradient = (1-gamma)*gradient + gamma*gradientnew;
    elseif strcmp(opts.methods,'recursiv')
        det1 = x(1)^2 - x(2:n+1)'*x(2:n+1);
        tmp1 = [x(1); -x(2:n+1)]/det1;
        det2 = x(n+2)^2 - x(n+3:2*n+2)'*x(n+3:2*n+2);
        tmp2 = [x(n+2); -x(n+3:2*n+2)]/det2;
        invx = [tmp1; tmp2];
        tmpgradient = data.train_miss(batch,:)'*(data.train_miss(batch,:)*x(2:n+1)-data.train_label(batch))/batchsize;
        tmpgradient = 2*tmpgradient./(1+tmpgradient.^2).^2;
        tmpgradient = [0; tmpgradient; zeros(size(tmpgradient,1)+1,1) ];
        lambdagradient1 = zeros(2*n+2,1);
        lambdagradient1(1) = lambda1;
        lambdagradient1(n+2) = lambda2;
        gradientnew = tmpgradient +  lambdagradient1;
        if i ==  1
            gradientbar = gradientnew;
        end
        gradientbar = gradientnew + (1 - gamma)* (gradientbar - gradientnew);
        gradient = gradientbar - muk * (invx + tmpgradient);
    end

    tmp{1} = 2*[x(1:n+1)*(x(1:n+1)'*gradient(1:n+1)) ; x(n+2:2*n+2)*(x(n+2:2*n+2)'*gradient(n+2:2*n+2))]  - [ det1*[gradient(1);-gradient(2:n+1)];det2*[gradient(n+2);-gradient(n+3:end) ]];
    lambda = - At{1}'*tmp{1};

    tmpschur1 = At{1}(1:n+1,:)'*x(1:n+1);
    tmpschur2 = At{1}(n+2:2*n+2,:)'*x(n+2:2*n+2);
    schur = 2*(tmpschur1*tmpschur1' + tmpschur2*tmpschur2') + (det1)*At{1}(1:n+1,:)'*At{1}(1:n+1,:) + det2*At{1}(n+2:2*n+2,:)'*At{1}(n+2:2*n+2,:);

    if ~issparse(schur); schur = sparse(schur); end;
    L.matfct_options = 'spchol';
    [L.R,indef,L.perm] = chol(schur,'vector');
    L.Rt  = L.R';
    L.matdim = length(schur);
    coeff.mat11 = schur;
    coeff.mat22 = [];
    coeff.mat12 = [];
    [lambda_2,resnrm,solve_ok] = symqmr(coeff,lambda,L,[],[],10);
    xold = x;
    tmplambda = At{1}*(lambda_2);

    tmpg2 = gradient + tmplambda;
    tmpgradient = 2*[x(1:n+1)*(x(1:n+1)'*tmpg2(1:n+1)) ;x(n+2:2*n+2)*(x(n+2:2*n+2)'*tmpg2(n+2:2*n+2))] - [ det1*[tmpg2(1);-tmpg2(2:9)];det2*[tmpg2(10);-tmpg2(11:end) ]];
        tmperrornew = sqrt(tmpg2'*tmpgradient);
        
        errortmp = norm(tmperrornew,'fro')^2;
        if i == 1
            tmperror0 = sqrt(errortmp);
        end
        if batchnum~= 0 
        tmperrorsum = tmperrorsum + errortmp;
        tmperror = sqrt(tmperrorsum/batchnum)/tmperror0 ;
        else
        tmperrorsum = tmperrorsum + errortmp;
        if ~(opts.options == 2)
        tmperror = sqrt(tmperrorsum/numall)/tmperror0 ;
        else
        tmperror = sqrt(tmperrorsum)/tmperror0  ;
        end
        tmperrorsum  = 0;
         end

    if ~strcmp(opts.methods,'recursiv')
        lr = opts.lr(i);
        % tmpnorm = [x(1:n+1)'*tmpg2(1:n+1);  x(1)*tmpg2(2:n+1)+tmpg2(1)*x(2:n+1)  ; x(n+2:2*n+2)'*tmpg2(n+2:2*n+2);  x(n+2)*tmpg2(n+3:2*n+2)+tmpg2(n+2)*x(n+3:2*n+2)  ];%inprod


        % if errortmp<1e-7
        %     break;
        % end
        ftmp = data.train_miss*x(2:n+1)-data.train_label;
        ftmp = sum(ftmp.^2./(1+ftmp.^2),"all")/size(data.train_miss,1) + x(1)*lambda1 + x(n+2)*lambda2;%
        out.error = [out.error tmperror];
        out.fval = [out.fval ftmp];
        x = xold - lr*(tmpgradient)/norm(tmperrornew,'fro');
    else
        tmperrornew = sqrt(tmpg2'*tmpgradient);
        lr = opts.lr(i)/norm(tmperrornew,'fro');
        tmperrornew = [x(1:n+1)'*tmpg2(1:n+1);  x(1)*tmpg2(2:n+1)+tmpg2(1)*x(2:n+1)  ; x(n+2:2*n+2)'*tmpg2(n+2:2*n+2);  x(n+2)*tmpg2(n+3:2*n+2)+tmpg2(n+2)*x(n+3:2*n+2)  ];%inprod

        ftmp = data.train_miss*x(2:n+1)-data.train_label;
        ftmp = sum(ftmp.^2./(1+ftmp.^2),"all")/size(data.train_miss,1)+ x(1)*lambda1 + x(n+2)*lambda2;
        out.error = [out.error tmperror];
        out.fval = [out.fval ftmp];
        x = xold - lr*tmpgradient;
    end


    muk = max(muk/1.2,1e-12);

end
out.x =x;

end