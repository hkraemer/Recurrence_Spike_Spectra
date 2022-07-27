% Here we compute and plot the actual corrcoeff between input signal and 
% re-composed signal as a function of the regularization parameter.

clear, clc

% generate some data
dt = 0.1;  % sampling
N = 1000; % length of considered tau-RR
NN = 2000; % length of considered trajectory
transients = 5000;
e = 0.1; % recurrence threshold
e_chaos = 0.05; % threshold for chaotic dynamics
thres_meth = 'var';
sigma = 0.05; % noise level

% parameters and initial condition
a = 0.36;
b = 2;
c = 4;
x0=[.7, -1, 0.4];

[~,x] = Rossler(a,b,c,x0,NN,dt,transients,0);

Y1 = x';
Y1 = (Y1 - mean(Y1)) ./ std(Y1);

RP1 = rp(Y1,e,thres_meth);
tau_rr1 = tau_recurrence_rate(RP1);

s = tau_rr1(1:200)';

%% Compute inter spike spectra

lambdas = 0:0.0005:0.05; % the encountered regularization thresholds

% normalize time series
s = (s - mean(s)) ./ std(s);
s = s - min(s);
s = s ./ max(s);

N = length(s);
% get set of basis functions
Theta = generate_basis_functions(N)';

coeffs_stls = zeros(1,length(lambdas));
coeffs_lasso = zeros(1,length(lambdas));


for i = 1:length(lambdas)
    display(i)
    actual_lambda = lambdas(i);
    % STLS
    y_act = stls(s, Theta, actual_lambda);
    rr = corrcoef(regenerate_signal(Theta, y_act), s);
    coeffs_stls(i) = rr(2);

    % LASSO
    y_act = lasso(Theta, s, 'Lambda', actual_lambda, 'Alpha', 1);
    rr = corrcoef(regenerate_signal(Theta, y_act), s);
    coeffs_lasso(i) = rr(2);
end

%% Plot results
lw = 2;
fs = 20;
sw = 60;

ylims = [0.84 1];

figure
subplot(121)
scatter(lambdas,coeffs_lasso,sw,'filled')
title("LASSO")
ylabel("\rho_{Pearson}")
xlabel("regularization parameter \alpha")
set(gca,'LineWidth',lw,'FontSize',fs)
ylim(ylims)
grid on

subplot(122)
scatter(lambdas,coeffs_stls,sw,'filled')
title("STLS")
ylabel("\rho_{Pearson}")
xlabel("regularization parameter \alpha")
set(gca,'LineWidth',lw,'FontSize',fs)
ylim(ylims)
grid on
