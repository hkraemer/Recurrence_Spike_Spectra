%% We look at the paradigmatic Roessler system in 3 different setups,
% period-2, period-3 and chaos. Here we only integrate the system in order
% to save the trajectories for later use. Systematically, we add an
% increasing amount of additional noise to it

clear, clc

dt = 0.05;  % sampling
dt_c = 0.1; % sampling for chaotic system
N = 1000; % length of considered tau-RR
NN = 2000; % length of considered trajectory
transients = 5000;
e = 0.1; % recurrence threshold
e_chaos = 0.05; % threshold for chaotic dynamics
thres_meth = 'var';
num_levels = 51; % number of noise levels
sigmas = linspace(0,0.5,num_levels); % noise levels

% parameters and initial condition
b = 2;
c = 4;
x0=[.7, -1, 0.4];

% preallocation of result vectors
tau_rr_1_all = zeros(num_levels, N);
tau_rr_2_all = zeros(num_levels, N);
tau_rr_3_all = zeros(num_levels, N);


%% Compute trajectories and corrresponding tau-RRs

% Limit-2
a = 0.36;
[~,x] = Rossler(a,b,c,x0,NN,dt,transients,0);
Y1 = x';
Y1 = (Y1 - mean(Y1)) ./ std(Y1);

% Limit-3
a = 0.41;
[~,x] = Rossler(a,b,c,x0,NN,dt,transients,0);
Y2 = x';
Y2 = (Y2 - mean(Y2)) ./ std(Y2);

% Chaos
rng(2) % random seed for reproducibility
x0 = randn(1,3);
a = 0.428;
[~,x] = Rossler(a,b,c,x0,NN,dt_c,transients,0);
Y3 = x';
Y3 = (Y3 - mean(Y3)) ./ std(Y3);

% compute tau RRs for different noise levels
for i = 1:length(sigmas)

    display(i)

    sigma = sigmas(i); % pick the noise level

    rng(123) % random seed for reproducibility
    Y1_n = Y1 + sigma*randn(size(Y1,1),size(Y1,2));
    RP1 = rp(Y1_n, e, thres_meth);
    tau_rr1 = tau_recurrence_rate(RP1);
    tau_rr1_all(i,:) = tau_rr1(1:N)/max(tau_rr1(1:N));

    rng(123) % random seed for reproducibility
    Y2_n = Y2 + sigma*randn(size(Y2,1),size(Y2,2));
    RP2 = rp(Y2_n, e, thres_meth);
    tau_rr2_n = tau_recurrence_rate(RP2);
    tau_rr2_all(i,:) = tau_rr2_n(1:N)/max(tau_rr2_n(1:N));

    rng(123) % random seed for reproducibility
    Y3_n = Y3 + sigma*randn(size(Y3,1),size(Y3,2));
    RP3 = rp(Y3_n, e_chaos, thres_meth);
    tau_rr3_n = tau_recurrence_rate(RP3);
    tau_rr3_all(i,:) = tau_rr3_n(1:N)/max(tau_rr3_n(1:N));

end



%% Save data

save("./computed data/tau_rr_1_all.csv", "tau_rr1_all","-ascii")
save("./computed data/tau_rr_2_all.csv", "tau_rr2_all","-ascii")
save("./computed data/tau_rr_3_all.csv", "tau_rr3_all","-ascii")
