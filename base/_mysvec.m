%% This function needs to be modified to support Var_sdp

% % This file is modified on SDPT3 svec.c

% compute the vector v = svec(X): extract the upper triangular part of a symmetric matrix M
% usage:
%   x = svec(K, M)
%   x = svec(K, M, isspM)
% input
%   K: cone object
%   M: one single matrix or a MatCell of matrices
%   is_sparse: either a integer or a vector of length length(X), indicating whether each element are treated as sparse or dense  -1: auto detect, 0: dense, 1: sparse


function x = svec(K, M, isspx)
     
    if nargin < 3  % auto detect 
        isspx = ones(length(K),1); 
    elseif length(isspx) == 1
        isspx = isspx * ones(length(K),1); 
    else
        assert(length(isspx) == length(K), 'isspx must be either a integer or of the same length as K');
    end

    %% M is one single matrix
    if ~ isa(M, 'MatCell') && ~ iscell(M) 
        if isa(K, 'Cone')
            assert(length(K) == 1);
            cone = K{1};
        elseif isa(K, 'BasicCone')
            cone = K;
        else
            error('K must be either a Cone or a BasicCone object');
        end
        
        if strcmp(cone.type,'s')
            if length(cone.size) > 1 && ~issparse(M)
                x = mexsvec({'s', cone.size}, sparse(M), 1); 
            else
                x = mexsvec({'s', cone.size}, sparse(M)); 
            end
        else
            x = M;
        end
        return
    end

    %% M is a MatCell consisting of matrices
    assert(length(K) == length(M), 'K and M must have the same length');
    x = MatCell(length(M)); 
    for p=1:length(K)
        cone = K{p};
        n = sum(cone.size);  
        Mp = M{p};
        if ~ isa(Mp, 'MatCell') && ~ iscell(Mp) % Mp is one single matrix
            Mp = {Mp}; 
        else
            % Mp is a MatCell concatenating several matrices
        end
        m = length(Mp);
        
        if strcmp(cone.type, 's')
            n2 = sum(cone.size .* (cone.size+1) ) / 2; 
            if isspx(p) 
                x{p} = sparse(n2,m); 
            else 
                x{p} = zeros(n2,m);  
            end

            for k = 1: m
                if length(cone.size) > 1 && ~issparse(M{p,k})  
                    assert(issymmetric(M{p,k}), 'M{p,k} must be a symmetric matrix');
                    x{p}(:,k) = mexsvec({'s', cone.size}, sparse(Mp{k}), isspx(p)); 
                else
                    x{p}(:,k) = mexsvec({'s', cone.size}, Mp{k}, isspx(p)); 
                end   
            end         
        else
            if isspx(p) 
                x{p} = sparse(n,m); 
            else 
                x{p} = zeros(n,m);  
            end
            for k = 1:m
                x{p}(:,k) = Mp{k}; 
            end
        end
    end
       
end 