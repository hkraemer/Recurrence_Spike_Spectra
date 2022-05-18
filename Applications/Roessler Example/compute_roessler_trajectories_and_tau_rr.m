%% We look at the paradigmatic Roessler system in 3 different setups,
% period-2, period-3 and chaos. Here we only integrate the system in order
% to save the trajectories for later use

clear, clc

dt = 0.05;  % sampling
dt_c = 0.1; % sampling for chaotic system
N = 1000; % length of considered tau-RR
NN = 2000; % length of considered trajectory
transients = 5000;
e = 0.1; % recurrence threshold
e_chaos = 0.05; % threshold for chaotic dynamics
thres_meth = 'var';
sigma = 0.05; % noise level

% parameters and initial condition
b = 2;
c = 4;
x0=[.7, -1, 0.4];

%% Limit-2
display("Compute Limit-2")

a = 0.36;
[~,x] = Rossler(a,b,c,x0,NN,dt,transients,0);

Y1 = x';
Y1 = (Y1 - mean(Y1)) ./ std(Y1);

RP1 = rp(Y1,e,thres_meth);
tau_rr1 = tau_recurrence_rate(RP1);
tau_rr1 = tau_rr1(1:N)/max(tau_rr1(1:N));

% 
%% Limit-3
display("Compute Limit-3")
a = 0.41;
[~,x] = Rossler(a,b,c,x0,NN,dt,transients,0);

Y2 = x';
Y2 = (Y2 - mean(Y2)) ./ std(Y2);

RP2 = rp(Y2,e,thres_meth);
tau_rr2 = tau_recurrence_rate(RP2);
tau_rr2 = tau_rr2(1:N)/max(tau_rr2(1:N));


%% Chaos
display("Compute chaos")

rng(2)
x0 = randn(1,3);
a = 0.428;
[~,x] = Rossler(a,b,c,x0,NN,dt_c,transients,0);

Y3 = x';
Y3 = (Y3 - mean(Y3)) ./ std(Y3);

RP3 = rp(Y3,e_chaos,thres_meth);
tau_rr3 = tau_recurrence_rate(RP3);
tau_rr3 = tau_rr3(1:N)/max(tau_rr3(1:N));


%% Noisy time series
display("Compute Limit-2 noise")
Y1_n = Y1 + sigma*randn(size(Y1,1),size(Y1,2));
RP1_n = rp(Y1_n,e,thres_meth);
tau_rr1_n = tau_recurrence_rate(RP1_n);
tau_rr1_n = tau_rr1_n(1:N)/max(tau_rr1_n(1:N));


%%
display("Compute Limit-3 noise")
Y2_n = Y2 + sigma*randn(size(Y2,1),size(Y2,2));
RP2_n = rp(Y2_n,e,thres_meth);
tau_rr2_n = tau_recurrence_rate(RP2_n);
tau_rr2_n = tau_rr2_n(1:N)/max(tau_rr2_n(1:N));


%%
display("Compute chaos noise")
Y3_n = Y3 + sigma*randn(size(Y3,1),size(Y3,2));
RP3_n = rp(Y3_n,e_chaos,thres_meth);
tau_rr3_n = tau_recurrence_rate(RP3_n);
tau_rr3_n = tau_rr3_n(1:N)/max(tau_rr3_n(1:N));


%% Save data

save("./computed data/tau_rr_1.csv", "tau_rr1","-ascii")
save("./computed data/tau_rr_1_n.csv", "tau_rr1_n","-ascii")
save("./computed data/tau_rr_2.csv", "tau_rr2","-ascii")
save("./computed data/tau_rr_2_n.csv", "tau_rr2_n","-ascii")
save("./computed data/tau_rr_3.csv", "tau_rr3","-ascii")
save("./computed data/tau_rr_3_n.csv", "tau_rr3_n","-ascii")

save("./computed data/Y_1.csv", "Y1","-ascii")
save("./computed data/Y_1_n.csv", "Y1_n","-ascii")
save("./computed data/Y_2.csv", "Y2","-ascii")
save("./computed data/Y_2_n.csv", "Y2_n","-ascii")
save("./computed data/Y_3.csv", "Y3","-ascii")
save("./computed data/Y_3_n.csv", "Y3_n","-ascii")

RP1 = RP1(1:1000,1:1000);
RP2 = RP2(1:1000,1:1000);
RP3 = RP3(1:1000,1:1000);
RP1_n = RP1_n(1:1000,1:1000);
RP2_n = RP2_n(1:1000,1:1000);
RP3_n = RP3_n(1:1000,1:1000);

save("./computed data/RP_1.csv", "RP1","-ascii")
save("./computed data/RP_1_n.csv", "RP1_n","-ascii")
save("./computed data/RP_2.csv", "RP2","-ascii")
save("./computed data/RP_2_n.csv", "RP2_n","-ascii")
save("./computed data/RP_3.csv", "RP3","-ascii")
save("./computed data/RP_3_n.csv", "RP3_n","-ascii")

