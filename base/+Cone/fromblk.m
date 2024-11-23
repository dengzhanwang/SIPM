%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-22 19:17:10
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-22 19:19:09
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function K = fromblk(blk)
    %% transform SDPT3 blk to Cone 
    K = {};
    for p = 1: size(blk, 1)
        cone = BasicCone(blk{p, 1}, blk{p, 2});
        K = [K, {cone}];    
    end
end
