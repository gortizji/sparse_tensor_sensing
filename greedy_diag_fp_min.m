function [select] = greedy_diag_fp_min(U,alpha,L,lambda)

%% Initialize variables

if nargin < 4
    lambda = 0;
end

R = length(U);  % Number of domains
select = cell(R,1); % Output selected indices
T_mem = cell(R,1);
N = zeros(R,1);
feasible_domains = 1:R;
empty_domains = [];
L_t = 0;

for i = 1:R
    N(i) = size(U{i},1); % Number of rows of Ui
    select{i} = 1:N(i);
    
    % Normalize rows of Ui to unit norm
    U{i} = normalize_rows(U{i}); 
    
    % Lookup table with inner products to speed up algorithm
    T_mem{i} = U{i}'*U{i};
    
    % Total number of selected indices
    L_t = L_t + length(select{i});
end

%% Greedy selection
FP_cand = zeros(R,1);

while L_t > L
    idx_aux = zeros(R,1); % Candidate indices to remove from set

    for i = feasible_domains
        
        % Initialize with first iteration
        % Remove nth row from Ui
        U_aux = U{i};
        U_aux(1,:) = [];
            
        % Update FP candidates
        Ti_aux = U_aux'*U_aux;
        FP_cand(i) = norm(Ti_aux.*hprod(T_mem, i), 'fro')^2 - regularizer(lambda, select, N, i);
        idx_aux(i) = 1;

        for n = 2:size(U{i},1)
            % Remove nth row from Ui
            U_aux = U{i};
            U_aux(n,:) = [];
            
            % Update FP candidates
            Ti_aux = U_aux'*U_aux;
            FP_aux = norm(Ti_aux.*hprod(T_mem, i), 'fro')^2 - regularizer(lambda, select, N, i);
            
            if FP_aux <= FP_cand(i)
                FP_cand(i) = FP_aux;
                idx_aux(i) = n;
            end
        end
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

    U_ = U{d_best};
    U_(idx_best,:) = [];
    U{d_best} = U_;
    T_mem{d_best} = U_'*U_;
        
    if length(select{d_best}) <= alpha(d_best)
        empty_domains(end+1) = d_best;
        temp_idx = find(feasible_domains == d_best);
        feasible_domains(temp_idx) = [];
    end
    L_t = L_t - 1;
end


end

function T_circ = hprod(T,j)
    N = length(T);
    T_circ = ones(size(T{1}));
    for i = 1:N
        if i ~= j
            T_circ = T_circ .* T{i};
        end
    end
end

function [T] = num_measurements(S, N, j)
    if nargin < 3
        j = -1;
    end
    T = 1;
    for i = 1:length(S)
        if i == j
            T = T * (N(i) + 1 - numel(S{i}));
        else
            T = T * (N(i) - numel(S{i}));
        end
    end
end

function [R] = regularizer(lambda, S, N, j)
if nargin < 4
    j = -1;
end

if lambda == 0
    R = 0;
else
    R = lambda * log(num_measurements(S, N, j));
end

end

