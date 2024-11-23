%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-22 19:19:29
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-26 10:54:37
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function blk = toblk(K)
%% return SDPT3 type representation
if iscell(K)
    blk = cell(length(K), 2);
    blk(:, 1) = cellfun(@(x) x.type, K, 'UniformOutput', false);
    blk(:, 2) = cellfun(@(x) x.size, K, 'UniformOutput', false);
else % K is a single BasicCone
    blk = cell(1, 2);
    blk{1, 1} = K.type;
    blk{1, 2} = K.size;
end
end