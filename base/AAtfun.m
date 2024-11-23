%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-18 17:07:34
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-11 21:04:17
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function Y = AAtfun(At)
% AAtfun  Compute A*A' for a given matrix A
if iscell(At) 
    m = size(At{1}, 2);
    Y = sparse(m, m);
    for p = 1: length(At)
        Y = Y + At{p}' * At{p};
    end
    pertdiag = 1e-13*ones(m,1); 
    Y = Y + spdiags(pertdiag,0,m,m);
else % At is a single matrix
    Y = At' * At;
end