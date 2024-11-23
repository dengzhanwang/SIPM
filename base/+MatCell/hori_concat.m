%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 14:10:16
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 09:47:14
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = hori_concat(X, idx)
    % horizontal concatenation 
    % input X is a cell

    if nargin < 2
        idx = 1:length(X);
    end
    out = horzcat(X{idx});
end