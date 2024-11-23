%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 08:49:54
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

%% M = mysvec(K, xvec, isspM)
function M = mysmat(varargin) 
    K = varargin{1};  % K is either Cone or BasicCone
    xvec = varargin{2};
    if iscell(xvec)
        M = MatCell(length(K));
        M.data = smat(Cone.toblk(K), xvec.data, varargin{3:end});
    else % xvec is one single matrix
        M = smat(Cone.toblk(K), xvec, varargin{3:end});
    end

end