classdef ViterbiDetector < handle
    properties(SetAccess=private)
        alphabet;
        M;
        K;
        psi;
        psi_delay;
        L1;
        L2;
    end
    
    properties(Access=private)
        states;
        connections;
        received_samples;
    end
    
    methods
        function self = ViterbiDetector(alphabet, trellis_depth, L1, L2, psi, psi_delay)
            if size(alphabet, 1) ~= 0
                alphabet = transpose(alphabet);
            end
            self.alphabet = alphabet;
            self.M = length(alphabet);
            self.K = trellis_depth;
            self.psi = psi;
            self.psi_delay = psi_delay;
            self.L1 = L1;
            self.L2 = L2;
            
            % Precompute some data
            self.build_all_states();
            self.build_all_connections();
            self.build_received_samples();
        end
        
        function [detected_syms, final_path_metrics] = detect_symbols(self, rho, initial_path_metrics)
            path_metrics = sparse(self.M^(self.L1+self.L2), self.K+1);
            survivors = sparse(self.M^(self.L1+self.L2), self.K+1);
            
            path_metrics(:,1) = initial_path_metrics; % First state is the initial state k=-1
            
            for k=0:self.K-1
                next_path_metrics = zeros(length(self.states), 1);
                next_survivors = zeros(length(self.states), 1);
                for j=1:length(self.states)
                    best_i = 0;
                    best_path_metric = Inf;
                    for i=1:length(self.states)
                        previous_path_metric = full(path_metrics(i, (k-1)+2));
                        if previous_path_metric == Inf; continue; end
                        branch_metr_i_j = self.branch_metric(j, i, rho(k+1));
                        if branch_metr_i_j == Inf; continue; end
                        next_path_metric = previous_path_metric + branch_metr_i_j;
                        if next_path_metric < best_path_metric
                            best_i = i;
                            best_path_metric = next_path_metric;
                        end
                    end
                    next_path_metrics(j) = best_path_metric;
                    next_survivors(j) = best_i;
                end
                path_metrics(:, k+2) = next_path_metrics;
                survivors(:, k+2) = next_survivors;
            end
            
            % Detected state sequence
            final_path_metrics = full(path_metrics(:, self.K+1));
            [~, min_j] = min(final_path_metrics);
            state_seq = zeros(self.K+1, self.L1+self.L2);
            state_seq(self.K-1+2,:) = self.states(min_j,:);
            k = self.K-1;
            j = min_j;
            while k > -1
                previous_state_i = full(survivors(j, k+2));
                previous_state = self.states(previous_state_i,:);
                state_seq((k-1)+2,:) = previous_state;
                k = k-1;
                j = previous_state_i;
            end
            
            % Detected symbols
            detected_syms = zeros(self.K, 1);
            k = self.K-1;
            while k >= 0
                detected_syms(k+1) = state_seq(k+2, self.L1+1);
                k = k-1;
            end
        end
        
        function bm = branch_metric(self, j, i, rho_k)
            if isempty(self.connections(j, i))
                bm = Inf;
            else
                bm = abs(rho_k - self.received_samples(j, i))^2;
            end
        end
    end
    
    methods(Access=private)
        function build_all_states(self)
            L = self.L1 + self.L2;
            self.states = self.build_all_states_recursive_helper(L);
        end
        
        function s = build_all_states_recursive_helper(self, L)
            if L > 1
                s_ = all_states(self.alphabet, L-1);
                s = [];
                for i=1:length(self.alphabet)
                    for j=1:size(s_,1)
                        r = [s_(j,:), self.alphabet(i)];
                        s = [s; r];
                    end
                end
            elseif L == 1
                for i=1:length(self.alphabet)
                    s(i,1) = self.alphabet(i);
                end
            else
                error('L<1');
            end
        end
        
        function build_all_connections(self)
            L = self.L1+self.L2;
            is = [];
            js = [];
            for i=1:length(self.states)
                sigma_i = self.states(i,:);
                for j=1:length(self.states)
                    sigma_j = self.states(j,:);
                    if all(sigma_j(2:L) == sigma_i(1:L-1))
                        is = [is; i];
                        js = [js; j];
                    end
                end
            end
            self.connections = sparse(js, is, ones(length(is), 1));
        end
        
        function build_received_samples(self)
            [js, is] = find(self.connections);
            us = zeros(length(is), 1);
            for n=1:length(is)
                sigma_j = self.states(js(n),:);
                sigma_i = self.states(is(n),:);
                a = [sigma_j(1,:), sigma_i(1,length(sigma_i))];
                u = conv(a, self.psi);
                us(n) = u(self.psi_delay+self.L1+1);
            end
            self.received_samples = sparse(js, is, us);
        end
    end
end
