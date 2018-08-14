 %% SCHEME COMPARISON: DENSE CORE
% MSE accuracy of greedy_fp, random
clear all; close all; clc;

N = [50, 60, 70];
K = [10, 20, 15];
alpha = [2, 2, 2];
Lmax = 180; % Maximum size to simulate
Lmin = sum(K) + sum(alpha); % Minimum size to simulate

Lstep = 1; % Step of increase in L
Nexp = 100; % Number of random experiments
Nrand = 100; % Number of random draws

L = Lmin:Lstep:Lmax;

%% Method comparison

Ntest = length(L);

% Initialize memory
mse_g_fp = zeros(Ntest, Nexp);
mse_rand = zeros(Ntest, Nexp, Nrand);

samples_g_fp = zeros(Ntest, Nexp);
samples_rand = zeros(Ntest, Nexp, Nrand);

for j = 1:Nexp
    fprintf('Simulation #%d\n', j)
    disp('----------------------------')
    tic
    U = cell(R,1);
    for r = 1:R 
        U{r} = randn(N(r),K);
    end
    
    parfor i = 1:Ntest
        for r=1:R
            alpha = ones(R,1)+1;
            alpha(r) = K;
            select = greedy_kron_fp_min(U, L(i), alpha);
            mse_new = MSE_diag(U, select);
            if mse_new < mse_g_fp(i,j) || r == 1
                mse_g_fp(i,j) = mse_new;
                samples_g_fp(i,j) = number_samples(select);
            end
        end
        
        for k = 1:Nrand
            alpha = ones(R,1);
            alpha(randi(R)) = K;
            select = random_kron_sampling(N, alpha, L(i));
            mse_rand(i,j,k) = MSE_diag(U, select);
            samples_rand(i,j,k) = number_samples(select);
        end
    end
    toc
end



%% Plot against sensors
T=10;
plot_g_fp = mean(abs(mse_g_fp), 2);
plot_g_fp_max = max(abs(mse_g_fp), [], 2);
plot_g_fp_min = min(abs(mse_g_fp), [], 2);
best = plot_g_fp(end);

close all
figure(1)
plot(L, 10*log10(plot_g_fp/best), 'LineWidth',2)
hold all
for p = [10,90]
    rand_value = prctile(abs(mse_rand), 100-p, 3);
    plot_rand = mean(rand_value,2);
    plot_rand_max = max(rand_value, [], 2);
    plot_rand_min = min(rand_value, [], 2);
    plot(L, 10*log10(plot_rand/best), 'LineWidth',2)
end

legend('greedy-fp', 'rand-10','rand-90')
title('Diagonal core comparison')
xlabel('L (number of sensors)')
ylabel('MSE')

%% Plot against samples
T =	30;
temp = [vec(samples_rand(:, T, :)), vec(mse_rand(:, T, :))];
temp = datasample(temp, 0.1 * size(temp,1));
figure
scatter(100*temp(:,1)/prod(N), 10*log10(temp(:,2)/best),100, '.', 'MarkerEdgeAlpha', 0.4, 'MarkerFaceAlpha', 0.4)
hold all
scatter(samples_g_fp(:, T)*100/prod(N), 10*log10(mse_g_fp(:, T)/best), 100,'.', 'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 1)
xlabel('Compression')
ylabel('MSE')
xtickformat('percentage')
set(gca,'xscale','log')
axis([1e-2, 1e2, -10, 90])
set(gca,'xticklabel',arrayfun(@(x) sprintf('%.2f %%',x),get(gca,'xtick'),'un',0))
legend('rand',  'g-fp')


