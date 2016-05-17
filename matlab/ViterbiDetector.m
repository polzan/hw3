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
        
        path_metrics;
        survivors;
        full_trellis;
    end
    
    methods
        function self = ViterbiDetector(alphabet, trellis_depth, L1, L2, psi, psi_delay)
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
            
            self.reset();
        end
        
        function reset(self, initial_path_metrics)
            if nargin < 2
                initial_path_metrics = zeros(length(self.states), 1);
            end
            
            self.path_metrics = zeros(self.M^(self.L1+self.L2), self.K+1);
            self.survivors = zeros(self.M^(self.L1+self.L2), self.K+1);
            
            self.path_metrics(:,1) = initial_path_metrics; % First state is the initial state k=-1
            self.full_trellis = -1;
        end
        
        function detected_sym = one_iteration(self, rho_k)
            next_path_metrics_j = [];
            next_path_metrics_v = [];
            
            next_survivors_j = [];
            next_survivors_v = [];
            
            for j=1:length(self.states)
                best_i = 0;
                best_path_metric = Inf;
                [~, is] = find(self.connections(j,:));
                previous_path_metrics = self.path_metrics(is, self.full_trellis+2);
                received_samples = self.received_samples(j, is);
                branch_metrics = abs(rho_k .* ones(1, length(received_samples)) - received_samples).^2;
                for ii=1:length(is)
                    i = is(ii);
                    previous_path_metric = previous_path_metrics(ii);
                    received_sample = received_samples(ii);
                    if previous_path_metric == Inf; continue; end
                    branch_metr_i_j = branch_metrics(ii);
                    next_path_metric = previous_path_metric + branch_metr_i_j;
                    if next_path_metric < best_path_metric
                        best_i = i;
                        best_path_metric = next_path_metric;
                    end
                end
                next_path_metrics_j = [next_path_metrics_j; j];
                next_path_metrics_v = [next_path_metrics_v; best_path_metric];
                
                next_survivors_j = [next_survivors_j; j];
                next_survivors_v = [next_survivors_v; best_i];
            end
            self.path_metrics(next_path_metrics_j, (self.full_trellis+1)+2) = next_path_metrics_v;
            self.survivors(next_survivors_j, (self.full_trellis+1)+2) = next_survivors_v;
            self.full_trellis = self.full_trellis + 1;
            
            if self.full_trellis == self.K-1
                final_path_metrics = self.path_metrics(:, self.K+1);
                [~, min_j] = min(final_path_metrics);
                
                % Follow the survivors to k=0
                j = min_j;
                k = self.K-1;
                while k > 0
                    j = self.survivors(j, k+2);
                    k = k-1;
                end
                
                detected_sym = self.states(j, self.L1+1);
                
                % Clear oldest
                self.path_metrics = self.path_metrics(:, 2:self.K+1);
                self.survivors = self.survivors(:, 2:self.K+1);
                self.full_trellis = self.full_trellis - 1;
            else
                detected_sym = NaN;
            end
        end
        
        function detected_symbols = detect(self, rho)
            detected_symbols = zeros(length(rho), 1);
            for k=0:length(rho)-1
                % Every 5000 symbols shift back the path metrics
                if mod(k, 5000) == 0
                    min_pm = min(min(self.path_metrics));
                    self.path_metrics = self.path_metrics - min_pm;
                end
                detected_symbols(k+1) = self.one_iteration(rho(k+1));
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
                s_ = self.build_all_states_recursive_helper(L-1);
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
            self.connections = full(sparse(js, is, ones(length(is), 1)));
        end
        
        function build_received_samples(self)
            [js, is] = find(self.connections);
            u = zeros(length(is), 1);
            % Pad psi when -L1..L2 is larger
            padded_psi_delay = self.psi_delay;
            padded_psi = self.psi;
            if -self.L1 + padded_psi_delay < 0
                pad_amount = self.L1 - padded_psi_delay;
                padded_psi = [zeros(pad_amount, 1); padded_psi];
                padded_psi_delay = padded_psi_delay + pad_amount;
            end
            if self.L2 + padded_psi_delay > length(padded_psi) - 1
                pad_amount = self.L2 + padded_psi_delay - length(padded_psi) + 1;
                padded_psi = [padded_psi; zeros(pad_amount, 1)];
            end
            eta = padded_psi((-self.L1:self.L2) + padded_psi_delay + 1);
            for n=1:length(is)
                sigma_j = self.states(js(n),:);
                sigma_i = self.states(is(n),:);
                a = [sigma_j(1,:), sigma_i(1,length(sigma_i))];
                u(n) = a * eta;
            end
            self.received_samples = full(sparse(js, is, u));
        end
    end
end
