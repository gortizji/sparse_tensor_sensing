clear; close all; clc;
addpath tensorlab;
N = 50; K = 10; Ix = 50; Iy = 60; M = 100;

azimuth_angles = linspace(5,85,K); % azimuth angles in degrees
idx = randperm(K);  
azimuth_angles(idx) = azimuth_angles;

elevation_angles = linspace(5,85,K); % azimuth angles in degrees
idx = randperm(K);  

elevation_angles(idx) = azimuth_angles;
lambda = 0.25;                  % wavelength of signal carriers
dx = 0.1;                       % sensor distance

for r = 1:K
  zar = exp(2*pi/lambda*sin(azimuth_angles(r)*pi/180)*dx*1i);
  for i = 1:Ix, A(i,r) = zar^(i-1); end
  zer = exp(2*pi/lambda*sin(elevation_angles(r)*pi/180)*dx*1i);
  for i = 1:Iy, E(i,r) = zer^(i-1); end
end
C = ((rand(M, K))>.5)*2-1;

attenuation_vector = rand(1, K);



%% Simulation parameters
L = [15,25,35,80,210]; % Number of sensors to place
SNR = -50:10:20; % SNRs to simulate
e_sim = zeros(length(SNR), length(L));
Nsim = 200;  % Number of simulations
Nrand = 50;  % Number of random sampling repetitions

%% Compute system matrices

Psi = cell(length(L),1);
Psi_inv = cell(length(L),1);
ns = cell(length(L),1);
T = cell(length(L),1);
for l = 1:length(L)
    select = greedy_diag_fp_min({A,E,C},[2,2,K],L(l));
    L_A = select{1}; L_E = select{2}; L_C = select{3};
    Psi{l} = kr({A(L_A,:),  E(L_E,:), C(L_C,:)});
    Psi_inv{l} = pinv(Psi{l});
    T{l} = Psi_inv{l}*Psi_inv{l}';
    ns{l} = number_samples(select);
end

%% Run simulation

for n = 1:length(SNR)
    SNR(n)
    sigma = 10^(-SNR(n)/20);

    for l = 1:length(L)
        for m = 1:Nsim
            S = zeros(K,N);
            
            for r = 1:K
              s = double(rand(N,1)>0.5)*2 - 1;
              S(r,:) = attenuation_vector(r)*s;
            end
            Y = Psi{l} * S;
            W = sigma*(1/sqrt(2)*randn(size(Y))+1j/sqrt(2)*randn(size(Y)));
            Y = Y + W;
            S_hat = Psi_inv{l} * Y;
            e_sim(n,l) = e_sim(n,l) + mse(S,S_hat);
        end
        e_sim(n,l) = e_sim(n,l)/(Nsim * N);
    end
end

%% Run random sampling

e_rand = zeros(Nrand, length(L)-1);
ns_rand = zeros(Nrand,length(L)-1);
sigma = 10^(30/20);
for n = 1:Nrand
    n
    for l=1:(length(L)-1)
        select = random_kron_sampling([Ix,Iy,M],[1,1,1],L(l));
        L_A = select{1}; L_E = select{2}; L_C = select{3};
        Psi = kr({A(L_A,:),  E(L_E,:), C(L_C,:)});
        Psi_inv = pinv(Psi);
        ns_rand(n,l) = number_samples(select);

        for m = 1:Nsim
            S = zeros(K,N);
            
            for r = 1:K
              s = double(rand(N,1)>0.5)*2 - 1;
              S(r,:) = attenuation_vector(r)*s;
            end
            Y = Psi * S;
            W = sigma*(1/sqrt(2)*randn(size(Y))+1j/sqrt(2)*randn(size(Y)));
            Y = Y + W;
            S_hat = Psi_inv * Y;
            e_rand(n,l) = e_rand(n,l) + mse(S,S_hat);
        end
        e_rand(n,l) = e_rand(n,l)/(Nsim * N);
    end
end


%% Plot results
close all
figure(1)
for l = 1:length(L)
     semilogy(SNR, e_sim(:,l))
     hold on
end
legend('15','25','35','80')

%% MSE computation 
function p = mse(S,S_hat)
    p = norm(S-S_hat)^2;
end