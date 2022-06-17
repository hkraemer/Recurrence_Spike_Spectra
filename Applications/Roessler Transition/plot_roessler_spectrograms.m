% We plot the spectrograms from the Roessler with transition
clear, clc
% load the data
% reconstructed
ISS1 = load("./results/results_Roessler_N_1000_ISS_0_85.csv");
ISS2 = load("./results/results_Roessler_N_1000_ISS_0_9.csv");
ISS3 = load("./results/results_Roessler_N_1000_ISS_0_95.csv");
ISS4 = load("./results/results_Roessler_N_1000_ISS_0_99.csv");
% true
ISS1t = load("./results/results_Roessler_N_1000_ISS_0_85_true.csv");
ISS2t = load("./results/results_Roessler_N_1000_ISS_0_9_true.csv");
ISS3t = load("./results/results_Roessler_N_1000_ISS_0_95_true.csv");
ISS4t = load("./results/results_Roessler_N_1000_ISS_0_99_true.csv");
% FFT spectrograms
FFT1 = load("./results/results_Roessler_N_1000_FFT_tau_rr_recon.csv");
FFT2 = load("./results/results_Roessler_N_1000_FFT_tau_rr_true.csv");
FFT3 = load("./results/results_Roessler_N_1000_FFT_time_series.csv");

%%
as = 0.36:0.0001:0.43;
dt = 0.2;
% according sampling frequency
fs = 1/dt;
% size of the considered τ-RR
tau_window = 500;
w2 = tau_window/2;
% make a time vector for the Fourier analysis
tt = 0:dt:(tau_window-1)*dt;

% plotting params
% cax_max = 0.013;
Fs = 22;
lw2 = 1;

figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
surf(as,tt(1:w2),FFT1')
view(2)
shading flat
% xline(13.9,'r--','LineWidth',lw2)
% xline(17,'r--','LineWidth',lw2)
% xline(34,'r--','LineWidth',lw2)
% xline(47,'r--','LineWidth',lw2)
% xline(56,'r--','LineWidth',lw2)
% xline(61,'r--','LineWidth',lw2)
colormap(parula)
% caxis([0 cax_max])
title('ISS1')
ylabel('Period [\Delta t]')
ylim([tt(1) tt(w2)])
xlim([as(1) as(end)])
set(gca,'FontSize',Fs)
set(gca,'YScale','log')
grid on

subplot(212)
surf(as,tt(1:w2),FFT3')
view(2)
shading flat
colormap(parula)
% caxis([0 cax_max])
title('ISS1 true')
ylabel('Period [\Delta t]')
ylim([tt(1) tt(w2)])
xlim([as(1) as(end)])
set(gca,'FontSize',Fs)
set(gca,'YScale','log')
grid on


figure('Units','normalized','Position',[.01 .01 .99 .99])
subplot(211)
surf(as,tt(1:w2),ISS4')
view(2)
shading flat
% xline(13.9,'r--','LineWidth',lw2)
% xline(17,'r--','LineWidth',lw2)
% xline(34,'r--','LineWidth',lw2)
% xline(47,'r--','LineWidth',lw2)
% xline(56,'r--','LineWidth',lw2)
% xline(61,'r--','LineWidth',lw2)
colormap(parula)
caxis([0 cax_max])
title('ISS1')
ylabel('Period [\Delta t]')
ylim([tt(1) tt(w2)])
xlim([as(1) as(end)])
set(gca,'FontSize',Fs)
set(gca,'YScale','log')
grid on

subplot(212)
surf(as,tt(1:w2),ISS4t')
view(2)
shading flat
colormap(parula)
caxis([0 cax_max])
title('ISS1 true')
ylabel('Period [\Delta t]')
ylim([tt(1) tt(w2)])
xlim([as(1) as(end)])
set(gca,'FontSize',Fs)
set(gca,'YScale','log')
grid on