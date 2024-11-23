%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-12 09:50:29
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 09:51:21
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = apply_fn(obj, other)
    % apply_fn(obj, other)
    out = elementwiseOperation(obj, other, @(fn, X) fn(X));
end