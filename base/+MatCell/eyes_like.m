function Y = eyes_like(X)
    Y = cell(size(X));
    for i = 1:length(X)
        Y{i} = eye((size(X{i})));
    end
end