clear
close all
clc

% Load dataset
load movielens.mat

%% Build users and movies graphs
G_u = gsp_graph(W_users);
G_m = gsp_graph(W_movies);

G_u = gsp_compute_fourier_basis(G_u);
G_m = gsp_compute_fourier_basis(G_m);

U2 = G_u.U;
U1 = G_m.U;

%% Estimate ground truth with GRALS
k = 10;
lambda_l = 1; 
lambda_w = 0.1;
lambda_h = 0.1;

Lh = lambda_l * G_m.L + lambda_w * eye(size(G_m.L));
Lw = lambda_l * G_u.L + lambda_w * eye(size(G_u.L));

[M,~,~] = grals(Lh,Lw,Otraining,M,k,Otest);

%% Set GFT bandwidth
K1 = 20;
K2 = 20;

U1_tilde = U1(:,1:K1);
U2_tilde = U2(:,1:K2);

%% Compute active query sample
select = greedy_kron_fp_min({U1_tilde,U2_tilde},100,[5,5]);

L1 = select{1};
L2 = select{2};

%% Reconstruct ratings
M_hat = reconstruct_sample_with_model(M,U1_tilde, U2_tilde, L1, L2);

% Build unobserved nodes mask for evaluation
W = ones(G_u.N,G_m.N);
W(L2,L1) = 0;

% Compute RMSE estimate
rmse = sqrt(norm(W.*M.*Otest-W.*M_hat.*Otest,'fro')^2/sum(Otest(:).*W(:)))
