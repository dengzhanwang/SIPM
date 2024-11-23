%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 20:55:15
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-11 20:55:21
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function C = rdivide(A, B)
    C = elementwiseOperation(A, B, @rdivide);
end