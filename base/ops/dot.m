
function C = dot(A, B)
    if iscell(A) && iscell(B)
        C = full(sum(cellfun(@(x, y) full(sum(sum(x .* y))), A, B)));
    else
        C = full(sum(sum(A .* B)));
    end

end