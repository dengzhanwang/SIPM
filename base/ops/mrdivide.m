%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 20:55:01
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-11 20:55:08
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function C = mrdivide(A, B)
    C = elementwiseOperation(A, B, @mrdivide);
end