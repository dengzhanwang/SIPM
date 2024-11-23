function Y = randn_like(X)
    Y = cell(size(X));
    for i = 1:length(X)
        Y{i} = randn((size(X{i})));
    end
end