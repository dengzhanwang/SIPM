%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 20:53:47
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-11 20:54:34
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function C = plus(A, B)
    C = elementwiseOperation(A, B, @plus);
end