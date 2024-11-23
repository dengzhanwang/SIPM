function C = elementwiseOperation(A, B, op)
    if iscell(A) && iscell(B)
        C = cellfun(@(a, b) op(a, b), A, B, 'UniformOutput', false);
    elseif iscell(A) 
        C = cellfun(@(a) op(a, B), A, 'UniformOutput', false);
    elseif iscell(B)
        C = cellfun(@(b) op(A, b), B, 'UniformOutput', false);
    else
        C = op(A, B);
    end
end