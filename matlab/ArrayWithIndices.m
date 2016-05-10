classdef ArrayWithIndices < handle
    properties(Access=protected)
        data;
        offset;
        autoPadding;
    end
    
    methods
        function obj = ArrayWithIndices(a, offset)
            if nargin < 1
                a = [];
            end
            if nargin < 2
                offset = 0;
            end
            obj.data = a;
            obj.offset = offset;
            obj.autoPadding = false;
        end
        
        function setAutoPadding(obj, val)
            obj.autoPadding = logical(val);
        end
        
        function val = isAutoPadding(obj)
            val = obj.autoPadding;
        end
        
        function setData(obj, ind, val)
            map_i = obj.reverse_map_indices(ind);
            d = obj.data;
            try
                d(map_i) = val;
            catch e
                if strcmp(e.identifier, 'MATLAB:badsubscript')
                    obj.alignOffset(min(ind));
                    d = obj.data;
                    map_i = obj.reverse_map_indices(ind);
                    d(map_i) = val;
                else
                    rethrow(e);
                end
            end
            obj.data = d;
        end
        
        function val = getData(obj, ind)
            map_i = obj.reverse_map_indices(ind);
            d = obj.data;
            try
                val = d(map_i);
            catch e
                if strcmp(e.identifier, 'MATLAB:badsubscript') && obj.autoPadding
                    [known_k, known_k_ord] = intersect(ind, obj.getIndices());
                    known = d(obj.reverse_map_indices(known_k));
                    [~, unknown_k_ord] = setdiff(ind, obj.getIndices());
                    out = zeros(length(ind),1);
                    out(known_k_ord) = known;
                    out(unknown_k_ord) = zeros(length(unknown_k_ord), 1);
                    val = out;
                else
                    rethrow(e);
                end
            end
        end
        
        function alignOffset(obj, newoff)
            if newoff < obj.offset
                diff = obj.offset - newoff;
                obj.data = [zeros(diff, 1); obj.data];
                obj.offset = newoff;
            elseif newoff > obj.offset
                diff =  newoff - obj.offset;
                d = obj.data;
                if any(d(1:diff) ~= 0)
                    error('Cannot drop nonzero elements');
                end
                obj.data = d(diff+1:length(d));
                obj.offset = newoff;
            end
        end
        
        function obj_t = translate(obj, delta)
            obj_t = ArrayWithIndices(obj.data, obj.offset + delta);            
        end
        
        function o = getOffset(obj)
            o = obj.offset;
        end
        
        function k = map_indices(obj, i)
            k = i + 1 + obj.offset;
        end
        
        function i = reverse_map_indices(obj,k)
            i = k - obj.offset + 1;
        end
        
        function indices = getIndices(obj)
            indices = (1:length(obj.data)) - 1 + obj.offset;
        end
        
        function d = getAll(obj)
            d = obj.data;
        end
        
        function obj_f = flip(obj)
            obj_f = ArrayWithIndices(flip(obj.data), - obj.offset - length(obj.data) + 1);
        end
    end
end
