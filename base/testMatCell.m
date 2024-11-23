%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-09 14:34:50
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
% Test script for MatCell class

% Create MatCell objects from cell arrays

clear; 
addpath_sdpt3;
t1_total = 0;
t2_total = 0;
t3_total = 0;

m = 10000;
n = 10000;

x1 = rand(m, n);
x2 = rand(m, n);
mc1 = MatCell({x1});
mc2 = MatCell({x2});
c1 = {x1};
c2 = {x2};
% mc1 = MatCell(c1);
% mc2 = MatCell(c2);

profile on;


% warm up
for i = 1: 10

mc_sum = mc1 + mc2;

op = @plus;
data_sum = op(x1, x2);

op = @plus;
c_sum = ops(c1, '+', c2);


end


% clock
for i = 1: 50

t1 = tic;
mc_sum = mc1 + mc2;
t1_total = t1_total + toc(t1);


t2 = tic ;
data_sum = x1 + x2;
t2_total = t2_total + toc(t2);


t3 = tic ;
op = @plus;
c_sum = ops(c1, '+', c2);
t3_total = t3_total + toc(t3);


end


t1_total
t2_total
t3_total

profile off;

% profile viewer;


% Check arithmetic operations results
