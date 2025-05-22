clear
close all;
test_dir = './';
seed = 2024;
rng(seed);
addpath(genpath(test_dir));
fvaltable = zeros(2,2,7);
errortable = zeros(2,2,7);
traintable = zeros(2,7,5);
testtable = zeros(2,7,5);
pall = [1 2];
n0all = [200 500];

p = 2;
k = 2;
n0 = n0all(k);
[blk, data, At ] = multitask_read2('winequality-red.csv','winequality-white.csv','ENB2012_data.xlsx','AirQualityUCI.xlsx','abalone.data',n0,p);


opts.lambda1 = 1e-2;
opts.lambda2 = 1e-2;
opts.batchsize = n0;

opts.maxiter = 10000;
Amap = @(x) At'*x;
ATmap = @(y) At*y;
opts.Amap = Amap;
opts.ATmap = ATmap;
opts.n = 8;
opts.gamma = 0.1;
x0 = zeros(opts.n,p*5);
opts.Amap = @(X) AXmap(X, blk, At);
opts.ATmap = @(y) Atymap(y, blk, At);
b = data.b;

opts.numall =5*p;
opts.m = 5*p;
Omega0 = 1/opts.m*eye(opts.m,opts.m);
opts.Omega0 = Omega0;
opts.methods = 'sto';
opts.lr = @(i) min(max(3/i^(1/1.5),0.01),1);
opts.options = 1;


out1 = sto_ipm_multi_task3(blk,x0,data,At,b,opts);




function [AX, AXorg] = AXmap(X, K, At, Lchol)
AX = AXfun_sdpnal(K,At,X);
end

function Aty = Atymap(y, K, At, Lchol)
Aty = Atyfun_sdpnal(K, At, y);
end
