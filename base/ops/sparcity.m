function out = sparcity(obj)
    % sparsity of the whole cell 
    if length(obj) == 0
        out = 0;
        return
    else
        s = sum(cellfun(@numel, obj));
        if s == 0
            out = 0;
        else
            out = nnz(obj) / s;
        end
    end
end