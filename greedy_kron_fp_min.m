function [select] = greedy_kron_fp_min(U,L,alpha)  
%GREEDY_KRON_FP_MIN Computes a near-optimal vertex selection for a product graph
%   Greedy minimization of the frame potential of the GFT of U=kron(U1,U2)
%   U: cell array containing {U1,U2}
%   L: number of nodes to select in the factor graphs
%   alpha: slack variable for the constraints -> Li >= Ki + alpha(i)
    

%% Initialize variables
R = length(U);  % Number of domains
K = zeros(R,1); % Size of domains
N = zeros(R,1);
select = cell(R,1); % Output selected indices

if nargin < 3
    alpha = zeros(R,1);
end  

G = cell(R,1);
FP_mem = zeros(R,1);
feasible_domains = 1:R;
empty_domains = [];
L_t = 0;
for i = 1:R
    N(i) = size(U{i},1); % Number of rows of Ui
    K(i) = size(U{i},2);
    
    select{i} = 1:N(i);
    
    % Normalize rows of Ui to unit norm
    U{i} = normalize_rows(U{i}); 
    
    % Lookup table with inner products to speed up algorithm
    G{i} = U{i}*U{i}';
    G{i} = G{i} - eye(size(G{i}));
    G{i} = G{i}.^2;
    
    % State memory
    FP_mem(i) = sum(sum(triu(G{i})));
    
    % Total number of selected indices
    L_t = L_t + length(select{i});
end


%% Greedy selection

while L_t > L
    idx_aux = ones(R,1); % Candidate indices to remove from set
    G_aux = G;
    FP_aux = FP_mem;
    FP_cand = (prod(FP_mem)) * ones(R,1) ;
    %L_t
    for i = feasible_domains
        [~,j] = max(sum(G{i}));
        idx_aux(i) = j;
        
        % Remove jth entry from lookup table
        Gi = G{i};
        Gi(:,j) = [];
        Gi(j,:) = [];
        G_aux{i} = Gi;
        
        % Update FP candidates
        FP_aux(i) = sum(sum(triu(G_aux{i})));
        %FP_cand(i) = (FP_aux(i) * prod(FP_mem) / FP_mem(i)) - regularizer(lambda, select, N, i);
        FP_cand(i) = (FP_aux(i) * prod(FP_mem) / FP_mem(i));
    end
    
    
    % Find candidate domain with smallest FP
    temp = min(FP_cand(feasible_domains));
    d_best = find(FP_cand == temp);
    d_best = d_best(1);
    % Remove idx from that domain
    idx_best = idx_aux(d_best);
 
    select_best = select{d_best};
    select_best(idx_best) = [];
    select{d_best} = select_best;
    
    % Update lookup table
    G{d_best} = G_aux{d_best};
    
    FP_mem(d_best) = FP_aux(d_best);
    
    if length(select{d_best}) <= (K(d_best) + alpha(d_best))
        empty_domains(end+1) = d_best;
        temp_idx = (feasible_domains == d_best);
        feasible_domains(temp_idx) = [];
    end
    L_t = L_t - 1;
end

end

function [T] = num_measurements(S, N, j)
    if nargin < 2
        j = -1;
    end
    T = 1;
    for i = 1:length(S)
        if i == j
            T = T * (N(i) - numel(S{i}) + 1);
        else
            T = T * (N(i) - numel(S{i}));
        end
    end
end

function R = regularizer(lambda, S, N, j)
if nargin < 4
    j = -1;
end

if lambda == 0
    R = 0;
else
    R = lambda * log(num_measurements(S, N, j));
end

end

