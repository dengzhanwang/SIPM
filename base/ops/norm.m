function out = norm(A)
    out = sqrt(sum(cellfun(@(x) norm(x, 'fro')^2, A)));
end