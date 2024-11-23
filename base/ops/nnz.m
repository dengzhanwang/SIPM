function out = nnz(obj)
    % sum of non-zero elements in cell array
    if length(obj) == 0
        out = 0;
        return
    else
        out = sum(cellfun(@nnz, obj));
    end
end