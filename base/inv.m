function out = inv(varargin)
    % Check if the number of inputs is 2 and if they are of the desired types
    if nargin == 2 && isa(varargin{1}, 'MatCell') && isa(varargin{2}, 'Cone')
        % Your custom implementation for MatCell and Cone
        out = MatCellinv(varargin{1}, varargin{2});
    else
        % Delegate to MATLAB's built-in inv function
        out = builtin('inv', varargin{:});
    end
end

function result = MatCellinv(X, K)
    assert(length(X) == length(K), 'MatCellinv: length(X) ~= length(K)')
    result = MatCell(length(K));
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 's')
            result{p} = inv(X{p});
        elseif strcmp(cone.type, 'q')
            result{p} = soc_ops(X{p}, 'inv', cone.size);
        elseif strcmp(cone.type, 'l')
            result{p} = 1 ./ X{p};
        elseif strcmp(cone.type, 'u')
            result{p} = zeros(size(X{p}));
        end
    end

end
