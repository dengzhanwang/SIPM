%% This function needs to be modified to support Var_sdp

%% This file is modified on SDPT3 smat.c
% compute the symmetric matrix X = smat(v) 
% usage:
%   M = smat(K, xvec)
%   M = smat(K, xvec, isspM)
% input
%   K: cone object
%   xvec: one single vector or a MatCell of vectors
%   is_sparse: either a integer or a vector of length length(X), indicating whether each element are treated as sparse or dense             0 or missing: dense, 1: sparse


function M = smat(K, xvec, isspM)

    if nargin < 3
        isspM = zeros(length(K),1); 
    elseif length(isspM) == 1
        isspM = isspM * ones(length(K),1); 
    else
        assert(length(isspM) == length(K), 'isspM must be either a integer or of the same length as K');
    end

    %% xvec is one single vector
    if ~ isa(xvec, 'MatCell') && ~ iscell(xvec)
        assert(length(K) == 1);
        cone = K{1};
        if strcmp(cone.type,'s')
            M = mexsmat({'s', cone.size}, xvec, isspM);      
        else
            M = xvec; 
        end   
        return;
    end
    

    %% xvec is MatCell consisting of numbers of vectors
    assert(length(K) == length(xvec), 'K and xvec must have the same length');
    M = MatCell(length(K));
    for p=1:length(K)
        cone = K{p};
        if strcmp(cone.type,'s')
            M{p} = mexsmat({'s', cone.size}, xvec{p}, isspM(p));
        else
            M{p} = xvec{p}; 
        end   
    end
end
 
 
 