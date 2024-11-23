


%% In matlab, cell is more efficient than class.
%  However,  in C++ we recommend to use class to represent cone.


function K = Cone(cones)
    K = cones;
end


%% original class

% classdef  Cone 
%     properties  (Access = public)
%         data  % cell 
%         blk
%     end

%     methods 
%         function obj = Cone(varargin)
%             % Matcell(c) convert a cell to MatCell, where c is a cel containing BasicCones; 

%             if nargin == 1
%                 if iscell(varargin{1})
%                     assert (size(varargin{1}, 2) == 1 || size(varargin{1}, 1) == 1);
%                     % reshape cell to be [length, 1]
%                     obj.data = reshape(varargin{1}, [], 1);
%                 else
%                     error('Invalid input argument');
%                 end
%             else
%                 error('Invalid input argument');
%             end

%             obj.blk = Cone.toblk(obj);
%         end


% %% operator overload
%         function n = length(obj)
%             n = length(obj.data);
%         end

% %% indexing
%         function out = subsref(obj, s)
%             if strcmp(s(1).type, '{}')
%                 if length(s) == 1
%                     if isnumeric(s(1).subs{:}) && numel(s(1).subs{:}) > 1  % handle array indexing
%                         out = cell(1, numel(s(1).subs{:}));
%                         for i = 1:numel(s(1).subs{:})
%                             out{i} = obj.data{s(1).subs{i}};
%                         end
%                     else
%                         out = obj.data{s(1).subs{:}};
%                     end
%                 else
%                     out = subsref(obj.data{s(1).subs{:}}, s(2:end));
%                 end
%             else
%                 out = builtin('subsref', obj, s);
%             end
%         end

%         function obj = subsasgn(obj, s, val)
%             if strcmp(s(1).type, '{}')
%                 if length(s) == 1
%                     if isnumeric(s(1).subs{:}) && numel(s(1).subs{:}) > 1  % handle array indexing
%                         assert(iscell(val), 'Right hand side of assignment must be a cell array when indexing with an array');
%                         assert(numel(s(1).subs{:}) == numel(val), 'Number of indices must match number of values to assign');
%                         for i = 1:numel(s(1).subs{:})
%                             obj.data{s(1).subs{i}} = val{i};
%                         end
%                     else
%                         obj.data{s(1).subs{:}} = val;
%                     end
%                 else
%                     obj.data{s(1).subs{:}} = subsasgn(obj.data{s(1).subs{:}}, s(2:end), val);
%                 end
%             else
%                 obj = builtin('subsasgn', obj, s, val);
%             end
%         end

%         function disp(obj)
%             print_tol = 10;
%             fprintf(' Cone object of %d blocks\n', length(obj));
%             for p = 1: length(obj)
%                 cone = obj.data{p};
%                 fprintf("\tblock %i: '%s' ", p, cone.type);
%                 if length(cone.size) <= print_tol  % print all if not too many
%                     fprintf("[ ");
%                     for j = 1: length(cone.size)
%                         fprintf("%i ", cone.size(j));
%                     end
%                     fprintf("]\n");
%                 else    % print first and last few if too many
%                     fprintf("[ ");
%                     n_head = ceil(print_tol / 2); n_tail = print_tol - n_head;
%                     for j = 1: n_head
%                         fprintf("%i ", cone.size(j));
%                     end
%                     fprintf("... ");
%                     for j = length(cone.size) - n_tail + 1: length(cone.size)
%                         fprintf("%i ", cone.size(j));
%                     end
%                     fprintf("]\n");
%                     fprintf("\t\t(Cartesian product of %i cones with sizes ranging from %i to %i)\n", ...
%                         numel(cone.size), min(cone.size), max(cone.size));
%                 end  
%             end      
%         end

%         function out = size(obj)
%         %% concate size of all cones
%             out = cellfun(@(x) sum(x.size), obj.data);
%         end
            
%     end

%     methods (Static)
%         function obj = fromblk(blk)
%             %% transform SDPT3 blk to Cone 
%             K = {};
%             for p = 1: size(blk, 1)
%                 cone = BasicCone(blk{p, 1}, blk{p, 2});
%                 K = [K, {cone}];    
%             end
%             obj = Cone(K);
%         end

%         function out = toblk(obj)
%         %% return SDPT3 type representation
%             out = cell(length(obj), 2);
%             out(:, 1) = cellfun(@(x) x.type, obj.data, 'UniformOutput', false);
%             out(:, 2) = cellfun(@(x) x.size, obj.data, 'UniformOutput', false);
%         end

%     end 


%     % from sedumi blk to cone

%     % 
% end