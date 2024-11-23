%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 20:56:08
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 15:49:24
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function C = transpose(A)
    if iscell(A)
        C = cell(size(A));
        for i = 1: numel(A)
            C{i} = A{i}.' ;
        end
    else
        C = A.';
    end
end