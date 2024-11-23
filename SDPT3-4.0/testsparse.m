clear
clc
m = 5000;
A = sparse(m,m);
for i = 1:m
    A(i,i) = i;
end
B = randn(m,m);

tic 
Z1 = A*B;
% Z2 = randn(m,m)*Z1;
toc

tic
Z11 = full(A)*B;
% Z21 = Z11*randn(m,m);
toc