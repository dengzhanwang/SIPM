%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-26 10:55:26
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
% BasicCone dose not have to be a cone, it can also be some other set of constraints or indicating penalty in objective 
%% some common "cone" are listed below

% (1) mathematically well-defined cone
% free/unbounded:         type = 'u'; size = [n]; params = [];
% zeros:                  type = 'zero'; size = [n]; params = [];
% nonnegative orthant:    type = 'l'; size = [n]; params = [];
% nonpositive orthant:    type = 'lneg'; size = [n]; params = [];
% quadratic cone:         type = 'q'; size = [q1, q2, ..., qn]; params = [];  this cone is the cartesian product of n quadratic cones
% rotated quadratic cone: type = 'r'; size = [q1, q2, ..., qn]; params = []; this cone is the cartesian product of n rotated quadratic cones
% semidefinite cone:      type = 's'; size = [s1, s2, ..., sn]; params = []; this cone is the cartesian product of n semidefinite cones
% (2) other "cone" (not a cone, but a set of constraints)
% linear:                 type = 'box'; size = [n], params = [lb, ub]; where lu, ub are of size [n, 1]. This "cone" is not a cone, but a linear constraint lb <= x <= ub
% (3) other "cone" (not a cone, but indicating penalty in objective)
% lasso penalty:          type = 'lasso'; size = [n]; params = [lamda]; where lammca is a scalar. This "cone" is not a cone, but indicating a lasso penalty lamda * sum(abs(x)) in objective
% (4) user costomized cone


%% In matlab, struct is more efficient than class, so we use struct instead of class.
%  However,  in C++ we recommend to use class.
function obj = BasicCone(type, size, params)
    if nargin < 3
        params = [];
    end
    obj = struct();
    obj.type = type;
    obj.size = size;
    obj.params = params;
    %% return SDPT3 type representation
    % obj.blk = {obj.type, obj.size};
    
end




% classdef  BasicCone
%     properties
%         type
%         size
%         params
%     end

%     methods
%         function obj = BasicCone(type, size, params)
%             if nargin < 3
%                 params = [];
%             end
%             obj.type = type;
%             obj.size = size;
%             obj.params = params;
%         end

%         function out = blk(obj)
%             %% return SDPT3 type representation
%             out = {obj.type, obj.size};
%         end
        
%     end

% end

