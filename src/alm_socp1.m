function [out] = alm_socp1(x0,At,data,opts)

maxiter = opts.max_iter;  
rho = opts.rho;        
eta = 5e-1;       
alpha = opts.alpha;      
n = opts.n;          
lambda1 = opts.lambda1;    
lambda2 = opts.lambda2;    
batchsize = opts.batchsize;
out.error = [];
x = x0;
lambda = zeros(size(At{1},2),1);
g = zeros(size(x));


out.fval = [];

for k = 1:maxiter

    grad_x  = g + At{1}*lambda + rho *  At{1}*(At{1}'*x ) ;
    x = x - eta * grad_x;
    x(1:n+1) = projection_to_soc(x(1:n+1));
    x(n+2:2*n+2) = projection_to_soc(x(n+2:2*n+2));
    lambda = lambda + rho * (At{1}'*x);
    tmpgradient = data.train_miss'*(data.train_miss*x(2:n+1)-data.train_label);
    tmpgradient = 2*tmpgradient./(1+tmpgradient.^2).^2;
    if size(data.train_miss,1) == 1e4
    tmpgradient = tmpgradient*20;% ./(1+tmpgradient.^2).^2
    end
    tmpgradient = [0; tmpgradient; zeros(size(tmpgradient,1)+1,1) ];
    lambdagradient1 = zeros(2*n+2,1);

    gradient = tmpgradient +  lambdagradient1;
    g = alpha * gradient + (1 - alpha) * g ;
    ftmp = objective(x,data,batchsize,lambda1,lambda2,n);
    out.fval = [out.fval ftmp];
end

out.x = x;

end

function ftmp = objective(x,data,batchsize,lambda1,lambda2,n)
ftmp = data.train_miss*x(2:n+1)-data.train_label;
norm(ftmp,'fro');
 ftmp = sum(ftmp.^2./(1+ftmp.^2),"all")/+ x(1)*lambda1 + x(n+2)*lambda2;
end


function x_proj = projection_to_soc(x)
if norm(x(2:end)) > x(1)
    tmp = abs((x(1) + norm(x(2:end)))/2/norm(x(2:end)));
    x_proj = tmp*[norm(x(2:end));x(2:end)] ;
else
    x_proj = x;
end


end
