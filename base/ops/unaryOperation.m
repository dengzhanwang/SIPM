function C = unaryOperation(A, op)
    if iscell(A)
        C = cellfun(@(a) op(a), A, 'UniformOutput', false);
    else
        C = op(A);
    end
end