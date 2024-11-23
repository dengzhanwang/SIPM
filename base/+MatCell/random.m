%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 21:06:57
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 15:54:54
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function Y = random(X)
    Y = cell(size(X));
    for i = 1:length(X)
        Y{i} = randn((size(X{i})));
    end
end