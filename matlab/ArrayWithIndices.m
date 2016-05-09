classdef ArrayWithIndices
    properties(Access=protected)
        data;
    end
    
    properties
        offset;
        autoPadding;
    end
    
    properties(Dependent=true)
        indices;
        all;
    end
    
    methods
        function obj = ArrayWithIndices(a)
            if nargin < 1
                a = [];
            end
            obj.data = a;
            obj.offset = 0;
            obj.autoPadding = false;
        end
        
        function obj = subsasgn(obj, s, varargin)
            switch s(1).type
                case '.'
                    % Use built-in for any other expression
                    obj = builtin('subsasgn',obj,s,varargin{:});
                case '()'
                    if length(s) == 1
                        % Implement obj(indices) = varargin{:};
                        d = obj.data;
                        map_i = obj.reverse_map_indices(s.subs{1});
                        d(map_i) = varargin{1};
                        obj.data = d;
                    else
                        % Use built-in for any other expression
                        obj = {builtin('subsasgn',obj,s,varargin)};
                    end
                case '{}'
                    % Use built-in for any other expression
                    obj = {builtin('subsasgn',obj,s,varargin)};
                otherwise
                    error('Not a valid indexing expression')
            end
        end
        
        function varargout = subsref(obj,s)
            switch s(1).type
                case '.'
                    % Use built-in for any other expression
                    varargout = {builtin('subsref',obj,s)};
                case '()'
                    if length(s) == 1
                        % Implement obj(indices)
                        d = obj.data;
                        map_i = obj.reverse_map_indices(s.subs{1});
                        try
                            varargout = {d(map_i)};
                        catch e
                            if strcmp(e.identifier, 'MATLAB:badsubscript')
                                [known_k, known_k_ord] = intersect(s.subs{1}, obj.indices);
                                known = d(obj.reverse_map_indices(known_k));
                                [~, unknown_k_ord] = setdiff(s.subs{1}, obj.indices);
                                out = [];
                                out(known_k_ord) = known;
                                out(unknown_k_ord) = zeros(length(unknown_k_ord), 1);
                                varargout = {out};
                            else
                                rethrow(e);
                            end
                        end
                    else
                        % Use built-in for any other expression
                        varargout = {builtin('subsref',obj,s)};
                    end
                case '{}'
                    % Use built-in for any other expression
                    varargout = {builtin('subsref',obj,s)};
                otherwise
                    error('Not a valid indexing expression')
            end
            
        end
        
        function k = map_indices(obj, i)
            k = i + 1 + obj.offset;
        end
        
        function i = reverse_map_indices(obj,k)
            i = k - obj.offset + 1;
        end
        
        function indices = get.indices(obj)
            indices = (1:length(obj.data)) - 1 + obj.offset;
        end
        
        function d = get.all(obj)
            d = obj.data;
        end
        
        function obj_f = flip(obj)
            obj_f = ArrayWithIndices(flip(obj.data));
            obj_f.offset = - obj.offset - length(obj.data) + 1;
        end
    end
end
