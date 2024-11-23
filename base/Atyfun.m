%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-11 21:04:02
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

%% comput Aty = A'*y

function Aty = Atyfun(K, At, y)

if iscell(At) 
    Aty = MatCell(length(K));
    for p = 1:length(K)
        cone = K{p};
        n = sum(cone.size);
        if strcmp(cone.type,'s')
            if (isempty(At{p}))
                Aty{p} = sparse(n,n);
            else
                Aty{p} = mexsmat({cone.type, cone.size}, At{p}*y);
            end
        else
            if (isempty(At{p}))
                Aty{p} = zeros(n,1);
            else
                Aty{p} = At{p}*y;
            end
        end
    end
else% At is a single matrix
    cone = K; % K is a single cone
    n = sum(cone.size);
    if strcmp(cone.type,'s')
        if (isempty(At))
            Aty = sparse(n,n);
        else
            Aty = mexsmat_sdpnal({cone.type, cone.size}, At*y);
        end
    else
        if (isempty(At))
            Aty = sparse(n,n);
        else
            Aty = At*y;
        end
    end
end
%%*********************************************************

