% MatCell class is used store data in blocks
% MatCell can be used as cell equiped with blockwise operators, e.g.,+ - .* and so on.

% usage example see testMatCell.m


%% Note that once the MatCell is created, its length(i.e., the number of blocks)is fixed.




classdef  MatCell 
    properties  (Access = public)
        data  
        % data: if len == 1, we store data as a matrix, if len > 1, we store data as a cell 
        name
        % name: if len == 1, we store name as a matrix, if len > 1, we store name as a cell of char.
        % for len > 1 case 
        % User is able to access the block by name. 
        % This is not implemented yet. This feature will be convenient for ading slack variables when handling conic programming in general form.
        len
        % len: the number of blocks
    end

    methods 
        function obj = MatCell(varargin)
            % There are several king of constuctor:
            % MatCell(n) return an object of length n, where n is integer, each block is empty;
            % Matcell(c) convert a cell to MatCell, where c is a cell; 
            % Matcell(c, name) convert a cell to Matcell and the name of each block. The input can be a single matrix if user want to create a MatCell of length 1.
            % copy constructor

            if nargin == 1
                if iscell(varargin{1})
                    % MatCell(c) convert a cell to MatCell, where c is a cell;
                    assert (size(varargin{1}, 2) == 1 || size(varargin{1}, 1) == 1);
                    % reshape cell to be [length, 1]
                    obj.len = length(varargin{1});
                    if obj.len == 1
                        obj.data = varargin{1}{1};
                        obj.name = '';
                    else
                        obj.data = reshape(varargin{1}, [], 1);
                        obj.name = cell(obj.len, 1);
                    end
                elseif isnumeric(varargin{1}) && isscalar(varargin{1})
                    % MatCell(n) return an object of length n, where n is integer, each block is empty;
                    obj.len = varargin{1};
                    if obj.len == 1
                        obj.data = [];
                        obj.name = '';
                    else
                        obj.data = cell(obj.len, 1);
                        obj.name = cell(obj.len, 1);
                    end
                else
                    error('Invalid input argument');
                end
            elseif nargin == 2
                if iscell(varargin{1}) && iscell(varargin{2})
                    % Matcell(c, name) convert a cell to Matcell and the name of each block. The two inputs are both cell.
                    assert(length(varargin{1}) == length(varargin{2}));
                    obj.len = length(varargin{1});
                    if obj.len == 1
                        obj.data = varargin{1}{1};
                        obj.name = varargin{2}{1};
                    else
                        obj.data = reshape(varargin{1}, [], 1);
                        obj.name = reshape(varargin{2}, [], 1);
                    end
                elseif isnumeric(varargin{1}) && ischar(varargin{2})
                    % Matcell(c, name) convert a cell to Matcell and the name of each block. The first input is a single matrix and the second input is a char.
                    obj.len = 1;
                    obj.data = varargin{1};
                    obj.name = varargin{2};
                else
                    error('Invalid input argument');
                end
            else
                error('Too many input arguments');
            end

            assert(obj.len > 0, 'empty MatCell not allowed')
        end

%% overload copy function as deep copy
        function new_obj = copy(obj)
            new_obj = MatCell(obj.data, obj.name);
        end


%% operator overload

%% length
        function n = length(obj)
            n = obj.len;
        end

%% indexing
        function out = subsref(obj, s)
            if strcmp(s(1).type, '{}')
                subs = s(1).subs{:};
                if length(s) == 1
                    if isempty(subs)  % handle empty array indexing
                        out = [];
                        return;
                    end
                    if numel(subs) > 1  % handle array indexing
                        assert(numel(subs) <= obj.len, 'Index exceeds matrix dimensions');
                        assert(all(subs) <= obj.len, 'Index exceeds matrix dimensions');
                        out = cell(1, numel(subs));
                        for i = 1:numel(subs)
                            out{i} = obj.data{s(1).subs{i}};
                        end 
                    else
                        if obj.len == 1
                            assert(isequal(subs, 1));
                            out = obj.data;
                        else
                            out = obj.data{subs};
                        end
                    end
                else
                    if obj.len == 1
                        assert(isequal(subs, 1));
                        out = subsref(obj.data, s(2:end));
                    else
                        out = subsref(obj.data{subs}, s(2:end));
                    end
                end
            elseif strcmp(s(1).type, '()')
                subs = s(1).subs{:};
                if obj.len == 1
                    assert(isequal(subs, obj.name));
                    out = obj.data;
                else
                    assert(ischar(subs), "name must be char");
                    idx = find(strcmp(obj.name,subs));
                    if isempty(idx)
                        error('Invalid block name');
                    end
                    if length(s) == 1
                        out = obj.data{idx};
                    else
                        out = subsref(obj.data{idx}, s(2:end));
                    end
                end
            else
                out = builtin('subsref', obj, s);
            end
        end

        function obj = subsasgn(obj, s, val)
            if strcmp(s(1).type, '{}')
                if length(s) == 1
                    if isnumeric(s(1).subs{:}) && numel(s(1).subs{:}) > 1  % handle array indexing
                        assert(iscell(val), 'Right hand side of assignment must be a cell array when indexing with an array');
                        assert(numel(s(1).subs{:}) == numel(val), 'Number of indices must match number of values to assign');
                        for i = 1:numel(s(1).subs{:})
                            obj.data{s(1).subs{i}} = val{i};
                        end
                    else
                        if obj.len == 1
                            assert(isequal(s(1).subs{:}, 1));
                            obj.data = val;
                        else
                            obj.data{s(1).subs{:}} = val;
                        end
                    end
                else
                    if obj.len == 1
                        assert(isequal(s(1).subs{:}, 1));
                        obj.data = subsasgn(obj.data, s(2:end), val);
                    else
                        obj.data{s(1).subs{:}} = subsasgn(obj.data{s(1).subs{:}}, s(2:end), val);
                    end
                end
            elseif strcmp(s(1).type, '()')
                assert(ischar(s(1).subs{:}), "name must be char");
                idx = find(strcmp(obj.name,s(1).subs{:}));
                if isempty(idx)
                    error('Invalid block name');
                end
                if length(s) == 1
                    if obj.len == 1
                        assert(isequal(s(1).subs{:}, obj.name));
                        obj.data = val;
                    else
                        obj.data{idx} = val;
                    end
                else
                    
                    obj.data{idx} = subsasgn(obj.data{idx}, s(2:end), val);
                end
            else
                obj = builtin('subsasgn', obj, s, val);
            end
        end


%% display
        function disp(obj)
            fprintf(' MatCell object of %d blocks\n', length(obj));
            if length(obj) == 1
                fprintf('\t{');
                if ~isempty(obj.name) && ~isempty(obj.name{k})
                    fprintf('%s: ', obj.name{k});
                end
                for j = 1: length(size(obj.data))
                    if j > 1
                        fprintf('x');
                    end
                    fprintf('%d', size(obj.data, j));
                end
                fprintf(' %s}\n', class(obj.data));
                return;
            end
            for k = 1: length(obj)
                fprintf('\t{');
                if ~isempty(obj.name) && ~isempty(obj.name{k})
                    fprintf('%s: ', obj.name{k});
                end
                for j = 1: length(size(obj.data{k}))
                    if j > 1
                        fprintf('x');
                    end
                    fprintf('%d', size(obj.data{k}, j));
                end
                fprintf(' %s}\n', class(obj.data{k}));
            end
        end


%% mathematical operations
        function out = elementwiseOperation(obj, other, operation)
            out = MatCell(length(obj));
            if isa(obj, 'MatCell') && isa(other, 'MatCell') 
                assert(length(obj) == length(other), 'MatCell operator: length not match')
                if length(obj) == 1
                    out.data = operation(obj.data, other.data);
                    out.name = obj.name;
                else
                    for k = 1: length(obj)
                        out.data{k} = operation(obj.data{k}, other.data{k});
                        out.name{k} = obj.name{k};
                    end
                end
            elseif isa(obj, 'MatCell') 
                if length(obj) == 1
                    out.data = operation(obj.data, other);
                    out.name = obj.name;
                else
                    for k = 1: length(obj)
                        out.data{k} = operation(obj.data{k}, other);
                        out.name{k} = obj.name{k};
                    end
                end
            elseif isa(other, 'MatCell')
                if length(other) == 1
                    out.data = operation(obj, other.data);
                    out.name = other.name;
                else
                    for k = 1: length(other)
                        out.data{k} = operation(obj, other.data{k});
                        out.name{k} = other.name{k};
                    end
                end
            else
                error('MatCell operator: type not supported')
            end
        end

        function out = plus(obj, other)
            % +
            out = elementwiseOperation(obj, other, @plus);
        end
        
        function out = minus(obj, other)
            % -
            out = elementwiseOperation(obj, other, @minus);
        end
        
        function out = times(obj, other)
            % .*
            out = elementwiseOperation(obj, other, @times);
        end
        
        function out = rdivide(obj, other)
            % ./
            out = elementwiseOperation(obj, other, @rdivide);
        end

        function out = mrdivide(obj, other)
            % /
            assert(isnumeric(other) && isscalar(other), 'MatCell operator /: type not supported');
            out = elementwiseOperation(obj, other, @mrdivide);
        end

        function out = apply_fn(obj, other)
            % apply_fn(obj, other)
            out = elementwiseOperation(obj, other, @(fn, X) fn(X));
        end

        function out = unaryOperation(obj, operation)
            out = MatCell(length(obj));
            if length(obj) == 1
                out.data = operation(obj.data);
            else
                for k = 1: length(obj)
                    out.data{k} = operation(obj.data{k});
                end
            end
        end

        function out = uminus(obj)
            % -obj
            out = unaryOperation(obj, @uminus);
        end

        function out = transpose(obj)
            % obj .'
            out = unaryOperation(obj, @transpose);
        end

        function out = ctranspose(obj)
            % obj '
            out = unaryOperation(obj, @ctranspose);
        end

        function out = full(obj)
            % full(obj)
            out = unaryOperation(obj, @full);
        end

        function out = sparse(obj)
            % sparse(obj)
            out = unaryOperation(obj, @sparse);
        end

        function out = spdiag(obj)
            % spdiag(obj)
            out = unaryOperation(obj, @spdiag);
        end

        function out = norm(obj)
            if length(obj) == 1
                out = norm(obj.data , 'fro');
            else
                out = sum(cellfun(@(x) norm(x, 'fro'), obj.data));
            end
        end

        function out = dot(obj, other)
            % dot(obj, other)

            % I an not sure if sum(sum(A .* B)); can be more efficient by writing mex function
            if length(obj) == 1
                out = sum(sum(obj.data .* other.data));
            else
                out = sum(cellfun(@(x, y) full(sum(sum(x .* y))), obj.data, other.data));
            end

        end


        function out = mtimes(obj, other)  
            % '*'
            % now support scalar * MatCell, MatCell * scalar
            assert(isnumeric(other) && isscalar(other) || isnumeric(obj) && isscalar(obj), 'MatCell operator *: type not supported');
            
            if isnumeric(other) && isscalar(other)
                out = MatCell(length(other));
                assert(isa(obj, 'MatCell'), 'MatCell operator *: type not supported');
                if length(obj) == 1
                    out.data = other * obj.data;
                else
                    out = MatCell(length(obj));
                    for k = 1: length(obj)
                        out.data{k} = other * obj.data{k};
                    end
                end
            elseif  isnumeric(obj) && isscalar(obj)
                assert(isa(other, 'MatCell'), 'MatCell operator *: type not supported');
                out = MatCell(length(other));
                if length(other) == 1
                    out.data = obj * other.data;
                else
                    for k = 1: length(other)
                        out.data{k} = obj * other.data{k};
                    end
                end
            else
                error('MatCell operator *: type not supported')
            end
        end

        function out = sparcity(obj)
            % sparsity of the whole MatCell 
            if length(obj) == 0
                out = 0;
                return
            elseif length(obj) == 1
                if numel(obj.data) == 0
                    out = 0;
                else
                    out = nnz(obj.data) / numel(obj.data);
                end
            else
                s = sum(cellfun(@numel, obj.data));
                if s == 0
                    out = 0;
                else
                    out = nnz(obj) / s;
                end
            end
        end

        function out = nnz(obj)
            % nnz(obj)
            if length(obj) == 0
                out = 0;
                return
            elseif length(obj) == 1
                out = nnz(obj.data);
            else
                out = sum(cellfun(@nnz, obj.data));
            end
        end


        function out = reduce_sum(obj)
            % compute the sum of all the blocks
            if length(obj) == 0
                out = [];
                return
            elseif length(obj) == 1
                out = obj.data;
                return
            end

            % old version implementation (using for loop)
            % firstBlockDims = []; % Initialize firstBlockDims to an empty array
            % for k = 1: length(obj)
            %     % Check if the current block has the same dimensions as the first block
            %     if isempty(firstBlockDims)
            %         firstBlockDims = size(obj.data{k});
            %     elseif ~isequal(firstBlockDims, size(obj.data{k}))
            %         error('All blocks must have the same dimensions for reduce_sum');
            %     end
        
            %     if k == 1
            %         out = obj.data{1};
            %     else
            %         out = out + obj.data{k};
            %     end
            % end

            % compute the sum of all the blocks
            firstBlockDims = size(obj.data{1}); % Get the dimensions of the first block
        
            % Check if all blocks have the same dimensions
            if ~all(cellfun(@(x) isequal(size(x), firstBlockDims), obj.data))
                error('All blocks must have the same dimensions for reduce_sum');
            end
        
            % Sum all the blocks together
            out = sum(cat(3, obj.data{:}), 3);
        end




        

%% concatenation
        function out = vert_concat(obj, idx)
            % vertical concatenation 
            % input can be block or cell

            if nargin == 1
                idx = 1:length(obj);
            end
            

            out = []; % output empty matrix if length(obj) == 0
            numRows = -1; % Initialize numRows to an invalid value
        
            for k = idx
                % Check if the current block has the same number of columns as the previous blocks
                if numRows == -1
                    numRows = size(obj.data{k}, 2);
                elseif numRows ~= size(obj.data{k}, 2)
                    error('All blocks must have the same number of columns for vertical concatenation');
                end
        
                if iscell(obj)
                    out = [out;
                            obj{k}];
                elseif isa(obj, 'MatCell')
                    out = [out;
                            obj.data{k}];
                end
            end
        end


        function out = hori_concat(obj, idx)
            % horizontal concatenation
            % input can be MatCell or cell
        
            out = []; % output empty matrix if length(obj) == 0
            numCols = -1; % Initialize numCols to an invalid value
            
            if nargin == 1
                idx = 1:length(obj);
            end
            
            for k = idx
                % Check if the current block has the same number of rows as the previous blocks
                if numCols == -1
                    numCols = size(obj.data{k}, 1);
                elseif numCols ~= size(obj.data{k}, 1)
                    error('All blocks must have the same number of rows for horizontal concatenation');
                end
        
                if iscell(obj)
                    out = [out, obj{k}];
                elseif isa(obj, 'MatCell')
                    out = [out, obj.data{k}];
                end
            end
        end


    end


    

%% constructor by splitting matrix
    methods (Static)
        function obj = zeros_like(A, dense_or_sparse)
            %% zeros_like - Create a MatCell object with the same size as the input matrix filled with zeros
            if nargin == 1
                dense_or_sparse = 'dense';
            end
            assert(isa(A, 'MatCell'), 'input must be MatCell')
            obj = MatCell(length(A));
            if length(A) == 1
                if strcmp(dense_or_sparse, 'dense')
                    obj.data = zeros(size(A.data));
                elseif strcmp(dense_or_sparse, 'sparse')
                    obj.data = sparse(size(A.data, 1), size(A.data, 2));
                end
            else
                for k = 1:length(A)
                    if strcmp(dense_or_sparse, 'dense')
                        obj.data{k} = zeros(size(A.data{k}));
                    elseif strcmp(dense_or_sparse, 'sparse')
                        obj.data{k} = sparse(size(A.data{k}, 1), size(A.data{k}, 2));
                    end
                end
            end
        end

        function obj = hori_split(A, cone_size)
            % hori_split - Split columns of the matrix A into blocks
            % Input:
            %   A - matrix to be split
            %   cone_size - array containing the number of columns for each block
            %
            % Output:
            %   obj - MatCell object containing the split blocks
        
            % Check if the sum of cone_size equals the number of columns in A
            assert(sum(cone_size) == size(A, 2), 'The sum of cone_size must equal the number of columns in A');
        
            % Initialize the MatCell object with the length of cone_size
            obj = MatCell(length(cone_size));
            if length(cone_size) == 1
                obj.data = A;
                return
            end
            % Split the columns of A into blocks and store them in the MatCell object
            col_start = 1;
            for k = 1:length(cone_size)
                col_end = col_start + cone_size(k) - 1;
                obj.data{k} = A(:, col_start:col_end);
                col_start = col_end + 1;
            end
        end

        function obj = vert_split(A, cone_size)
            % vert_split - Split rows of the matrix A into blocks
            % Input:
            %   A - matrix to be split
            %   cone_size - array containing the number of rows for each block
            %
            % Output:
            %   obj - MatCell object containing the split blocks
        
            % Check if the sum of cone_size equals the number of rows in A
            assert(sum(cone_size) == size(A, 1), 'The sum of cone_size must equal the number of rows in A');
        
            % Initialize the MatCell object with the length of cone_size
            obj = MatCell(length(cone_size));
            if length(cone_size) == 1
                obj.data = A;
                return
            end
            % Split the rows of A into blocks and store them in the MatCell object
            row_start = 1;
            for k = 1:length(cone_size)
                row_end = row_start + cone_size(k) - 1;
                obj.data{k} = A(row_start:row_end, :);
                row_start = row_end + 1;
            end
            
        end



    end
    
end





