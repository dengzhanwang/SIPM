%% 

%% Var_sdp is used to store the SDP variables and also data matrix
% In matrix form, it is the Cartesian product of numbers of suquare matrices.
% In vector form, it is the concatenation of the vectors.
% These two form can be converted to each other by the function mysmat and mysvec.

classdef Var_sdp < MatCell
    properties
        cone_size
    end

    methods
        function obj = Var_sdp(varargin)
            obj@MatCell(varargin); % Call the superclass constructor
            assert(all(cellfun(@(x) size(x, 1) == size(x, 2), obj.data)));
            obj.cone_size = cellfun(@(x) size(x, 1), obj.data);
        end

        function obj = smat(obj)
            % convert the vector form to the matrix form
            obj.data = mexsmat({'s', obj.cone_size}, obj.data);
        end

        function obj = mysvec(obj)
            % convert the matrix form to the vector form
            obj.data = mexsvec({'s', obj.cone_size}, obj.data); 
        end

        function obj = inv(obj)
            % compute inverse of each matrix
            obj.data = cellfun(@(x) inv(x), obj.data);
        end

    end

    methods (Static)
        function obj = zeros(cone_size)
            obj = Var_sdp(cellfun(@(x) zeros(x), cone_size, 'UniformOutput', false));
        end

        function obj = rand(cone_size)
            obj = Var_sdp(cellfun(@(x) rand_psd(x), cone_size, 'UniformOutput', false));
        end
    end
end