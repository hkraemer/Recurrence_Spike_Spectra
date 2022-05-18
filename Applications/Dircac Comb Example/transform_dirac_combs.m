%% Consider a Dirac-Comb and its FFT

clear, clc

rng(2)
N = 1000; % length of dirac comb
T = 100; % inter-spike period

dirac1 = zeros(1,N);
dirac2 = zeros(1,N);

dirac1(1:T:end) = 1;
dirac2(1:T:end) = rand(1,length(1:T:length(dirac2)));

% compute new spectra (0.9 rho-threshold, LASSO)
spectrum1 = inter_spike_spectrum(dirac1,'threshold',0.9);
spectrum2 = inter_spike_spectrum(dirac2,'threshold',0.9);
% compute periodograms
[p1,w1] = periodogram(dirac1,[],length(dirac1));
[p2,w2] = periodogram(dirac2,[],length(dirac1));

%% Plot DCs, Periodograms & Inter Spike Spectra

fs = 18;
lw2 = 2;

figure('Units','normalized','Position',[.01 .01 .99 .99]) 
subplot(2,3,1)
plot(1:N,dirac1,'LineWidth',3)
set(gca,'LineWidth',lw2)
set(gca,'FontSize',fs)
xlabel('time t [a.u.]')
ylabel('[a.u.]')
title('Perfect Dirac comb (DC)')
ylim([0 1])
grid on

subplot(2,3,4)
plot(1:N,dirac2,'LineWidth',3)
set(gca,'LineWidth',lw2)
set(gca,'FontSize',fs)
xlabel('time t [a.u.]')
ylabel('[a.u.]')
title('DC with random amplitudes')
ylim([0 1])
grid on

subplot(2,3,2)
plot((w1./(2*pi)),p1,'LineWidth',3)
set(gca,'LineWidth',lw2)
set(gca,'FontSize',fs)
xlabel('frequency [1/t]')
ylabel('Power/frequency')
xlim([0 0.2])
ylim([0 0.035])
title('Periodogram of perfect DC')
grid on

subplot(2,3,5)
plot((w2./(2*pi)),p2,'LineWidth',3)
set(gca,'LineWidth',lw2)
set(gca,'FontSize',fs)
xlabel('frequency [1/t]')
ylabel('Power/frequency')
xlim([0 0.2])
ylim([0 0.035])
title('Periodogram of DC with random amplitudes')
grid on

subplot(2,3,3)
plot(1./(1:N/2),spectrum1,'LineWidth',3)
set(gca,'LineWidth',lw2)
set(gca,'FontSize',fs)
xlabel('frequency [1/t]')
ylabel('rel. Power/frequency')
xlim([0 0.1])
ylim([0 1])
title('Inter Spike Spectrum of perfect DC')
grid on

subplot(2,3,6)
plot(1./(1:N/2),spectrum2,'LineWidth',3)
set(gca,'LineWidth',lw2)
set(gca,'FontSize',fs)
xlabel('frequency [1/t]')
ylabel('rel. Power/frequency')
xlim([0 0.1])
ylim([0 1])
title('Inter Spike Spectrum of randomized DC')
grid on

