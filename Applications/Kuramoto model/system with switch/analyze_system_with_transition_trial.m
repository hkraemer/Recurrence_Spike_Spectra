%% Visualize results for Kuramoto with transition

% Computations have been carried out on the cluster (windowed spectrum determination)
% in the script `analyze_system_with_transition_cl.jl`

clear, clc
% load the results:

% lambda_1 = 1e-4 (in Julia terms):
spectrum1 = load("./results/spectrum1_trial.csv");
% lambda_2 = 5e-5 (in Julia terms):
spectrum2 = load("./results/spectrum2_trial.csv");

%% Plotting

t1 = 1:size(spectrum1,1);
t2 = 1:size(spectrum1,2);

fs = 22;
lw2 = 1;

len = 90;

cax_max = 0.013;


figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(t1,t2,spectrum1')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
title("Evolutionary inter spike spectrogram of Kuramoto with switch")
ylabel("Period [log_{10} a.u.]")
set(gca,'Yscale','log')
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\lambda = 1e-4")
grid on

subplot(4,1,[3,4])
surf(t1,t2,spectrum2')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
ylabel("Period [log_{10} a.u.]")
xlabel("time [a.u.]")
set(gca,'Yscale','log')
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\lambda = 1e-5")
grid on



figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(t1,t2,spectrum1')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
title("Evolutionary inter spike spectrogram of Kuramoto with switch")
ylabel("Period [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\lambda = 1e-4")
grid on

subplot(4,1,[3,4])
surf(t1,t2,spectrum2')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
ylabel("Period [a.u.]")
xlabel("time [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\lambda = 1e-5")
grid on
