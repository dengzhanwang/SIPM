%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 14:10:16
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 09:48:16
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = vert_concat(X, idx)
    % vertical concatenation 
    % input X is a cell

    if nargin < 2
        idx = 1:length(X);
    end
    out = vertcat(X{idx});
end