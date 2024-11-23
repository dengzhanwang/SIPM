function out = sto_ipm_lin(blk,At,Csto,b,opts)
lr = opts.lr;
for i = 1:opts.maxiter
invxv = 1./x;
gradient = sum(Csto(:,:,batch),3)/batchsize + muk * invxv;
lambda = Amap*gradient.*X.^2;
schur = schurmat_sblk(blk,At,par,schur,p,X,X);
lambda = linesolve(schur,lambda);
xnew = x - lr*(X*Atmap*lambda*X)/norm(X*gradient + Atmap *lambda,'fro');
end


end