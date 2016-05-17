classdef MaxLogMapDetector < handle
    properties(SetAccess=private)
        alphabet;
        M;
        Kd;
        Kin;
        psi;
        psi_delay;
        L1;
        L2;
        state_length;
        state_count;
    end
    
    properties(Access=private)
        states;
        connections; % sigma_j <-> sigma_i
        received_samples; % u_k = f(sigma_j, sigma_i)
        
        last_forward_metrics;
        initial_backward_metrics;
    end
    
    methods
        function self = MaxLogMapDetector(alphabet, trellis_depth, max_input_length, L1, L2, psi, psi_delay)
            self.alphabet = alphabet;
            self.M = length(alphabet);
            self.Kd = trellis_depth;
            self.Kin = max_input_length;
            self.psi = psi;
            self.psi_delay = psi_delay;
            self.L1 = L1;
            self.L2 = L2;
            self.state_length = L1+L2;
            
            % Precompute some data
            self.build_all_states();
            self.state_count = length(self.states);
            self.build_all_connections();
            self.build_received_samples();
            
            self.reset();
        end
        
        function reset(self, initial_path_metrics, final_path_metrics)
            if nargin < 2
                initial_path_metrics = zeros(length(self.states), 1);
            end
            if nargin < 3
                final_path_metrics = zeros(length(self.states), 1);
            end
            self.initial_backward_metrics = final_path_metrics;
            self.last_forward_metrics = initial_path_metrics;
        end
        
        function detected_symbols = detect(self, rho)
            detected_symbols = zeros(length(rho), 1);
            detected_symbols(1:self.Kd-1) = NaN .* ones(self.Kd-1, 1);
            detected_symbols(length(rho)-self.Kd+2:length(rho)) = NaN .* ones(self.Kd-1, 1);
            next_sym = self.Kd-1;
            next_rho = 0;
            while next_rho + self.Kin < length(rho)
                rho_chunk = rho((0:self.Kin-1)+next_rho+1);
                [~, good_sym_chunk] = self.detect_block(rho_chunk);
                detected_symbols(next_sym+1:next_sym+length(good_sym_chunk)) = good_sym_chunk;
                next_rho = next_rho + self.Kin - 2*(self.Kd-1);
                next_sym = next_sym + length(good_sym_chunk);
            end
            
            if next_rho < length(rho)
                %warning('Using less symbols in the final chunk');
                rho_chunk = rho(next_rho+1:length(rho));
                [~, good_sym_chunk] = self.detect_block(rho_chunk);
                detected_symbols(next_sym+1:next_sym+length(good_sym_chunk)) = good_sym_chunk;
            end
            
            %             for k=0:length(rho)-1
            %                 % Every 5000 symbols shift up the path metrics
            %                 if mod(k, 5000) == 0
            %                     max_pm = full(max(max(self.path_metrics)));
            %                     self.path_metrics = self.path_metrics + max_pm;
            %                 end
            %                 detected_symbols(k+1) = self.one_iteration(rho(k+1));
            %             end
        end
        
        function [detected_symbols, detected_symbols_good] = detect_block(self, rho)
            fm = self.build_forward_metric(rho, self.last_forward_metrics);
            bm = self.build_backward_metric(rho, self.initial_backward_metrics);
            
            fm_good = fm(:, (self.Kd-1:length(rho)-1-self.Kd+1)+2);
            bm_good = bm(:, (self.Kd-1:length(rho)-1-self.Kd+1)+1);
            state_metric = fm_good + bm_good;
            
            self.last_forward_metrics = fm(:,length(rho)-1 - self.Kd+1 + 2);
            
            detected_symbols_good = zeros(length(rho)-2*(self.Kd-1), 1);
            for l=1:length(rho)-2*(self.Kd-1)
                [~, max_state_i] = max(state_metric(:, l));
                detected_symbols_good(l) = self.states(max_state_i, self.L1+1);
            end
            detected_symbols = [ ...
                NaN .* ones(self.Kd-1, 1); ...
                detected_symbols_good; ...
                NaN .* ones(self.Kd-1, 1) ...
                ];
        end
    end
    
    methods(Access=private)
        %         function update_ch_transition_metrics(self, rho_k, k)
        %             [js, is, uks] = find(self.received_samples);
        %             chtm = zeros(length(js), 1);
        %             for i=1:length(js)
        %                 chtm(i) = -abs(rho_k - uks(i))^2;
        %             end
        %             self.ch_transition_metrics();
        %         end
        
        function fm = build_forward_metric(self, rho, initial_forward_metrics)
            Kin = length(rho);
            fm = sparse(self.state_count, Kin+1); % k=-1..Kin-1
            fm(:,1) = initial_forward_metrics; % First state is the initial state k=-1
            for k=0:Kin-1 - self.Kd+1 % Skip the last Kd-1 -> not used
                next_path_metrics_j = [];
                next_path_metrics_v = [];
                for j=1:self.state_count
                    best_path_metric = -Inf;
                    [~, is] = find(self.connections(j,:));
                    previous_path_metrics = full(fm(is, (k-1)+2));
                    received_samples = full(self.received_samples(j, is));
                    branch_metrics = -abs(rho(k+1) .* ones(1, length(received_samples)) - received_samples).^2;
                    for ii=1:length(is)
                        i = is(ii);
                        previous_path_metric = previous_path_metrics(ii);
                        received_sample = received_samples(ii);
                        if previous_path_metric == -Inf; continue; end
                        branch_metr_i_j = branch_metrics(ii); % <- optimize by precomputing c_k(j|i) ?
                        next_path_metric = previous_path_metric + branch_metr_i_j;
                        if next_path_metric > best_path_metric
                            best_path_metric = next_path_metric;
                        end
                    end
                    next_path_metrics_j = [next_path_metrics_j; j];
                    next_path_metrics_v = [next_path_metrics_v; best_path_metric];
                end
                fm(next_path_metrics_j, k+2) = next_path_metrics_v;
            end
        end
        
        function bm = build_backward_metric(self, rho, initial_backward_metrics)
            Kin = length(rho);
            bm = sparse(self.state_count, Kin+1); % k=0..Kin
            bm(:,Kin+1) = initial_backward_metrics; % First state is the last state k=Kin
            k = Kin-1;
            while k >= 0 + self.Kd - 1 % First useful metrics at k=Kd-1
                next_path_metrics_i = [];
                next_path_metrics_v = [];
                for i=1:self.state_count
                    best_path_metric = -Inf;
                    [js, ~] = find(self.connections(:,i));
                    previous_path_metrics = full(bm(js, (k+1)+1));
                    next_received_samples = full(self.received_samples(js, i));
                    if k == Kin-1
                        branch_metrics = zeros(length(next_received_samples), 1);
                    else
                        branch_metrics = -abs(rho((k+1)+1) .* ones(length(next_received_samples),1) - next_received_samples).^2;
                    end
                    for jj=1:length(js)
                        j = js(jj);
                        previous_path_metric = previous_path_metrics(jj);
                        %received_sample = next_received_samples(jj);
                        if previous_path_metric == -Inf; continue; end
                        branch_metr_i_j = branch_metrics(jj);
                        next_path_metric = previous_path_metric + branch_metr_i_j;
                        if next_path_metric > best_path_metric
                            best_path_metric = next_path_metric;
                        end
                    end
                    next_path_metrics_i = [next_path_metrics_i; i];
                    next_path_metrics_v = [next_path_metrics_v; best_path_metric];
                end
                bm(next_path_metrics_i, k+1) = next_path_metrics_v;
                k = k-1;
            end
        end
        
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
            self.connections = sparse(js, is, ones(length(is), 1));
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
            self.received_samples = sparse(js, is, u);
        end
    end
end
