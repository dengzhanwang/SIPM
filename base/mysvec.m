%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 16:57:37
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 08:50:05
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

%% x = mysvec(K, M, isspM)
function x = mysvec(varargin)
    K = varargin{1};
    M = varargin{2};
    if iscell(M)
        x = MatCell(length(K));
        x.data = mexsvec(Cone.toblk(K), sparse(M.data), varargin{3:end});
    else % M is one single matrix
        x = mexsvec(Cone.toblk(K), sparse(M), varargin{3:end});
    end

end