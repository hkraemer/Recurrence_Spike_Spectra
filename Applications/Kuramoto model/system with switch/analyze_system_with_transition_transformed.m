%% Visualize results for Kuramoto with transition

% Computations have been carried out on the cluster (windowed spectrum determination)
% in the script `analyze_system_with_transition_cl.jl`

clear, clc
% load the results for the different rho-thresholds:

spectrum1 = load("./results/results_kuramoto_switchspectrum_0_8.csv");
spectrum2 = load("./results/results_kuramoto_switchspectrum_0_85.csv");
spectrum3 = load("./results/results_kuramoto_switchspectrum_0_9.csv");
spectrum4 = load("./results/results_kuramoto_switchspectrum_0_95.csv");

spectrum1t = load("./results/results_kuramoto_switchspectrum_0_8_t.csv");
spectrum2t = load("./results/results_kuramoto_switchspectrum_0_85_t.csv");
spectrum3t = load("./results/results_kuramoto_switchspectrum_0_9_t.csv");
spectrum4t = load("./results/results_kuramoto_switchspectrum_0_95_t.csv");

rhos1 = load("./results/results_kuramoto_switchrhos1.csv");
rhos2 = load("./results/results_kuramoto_switchrhos2.csv");
rhos3 = load("./results/results_kuramoto_switchrhos3.csv");
rhos4 = load("./results/results_kuramoto_switchrhos4.csv");

rhos1t = load("./results/results_kuramoto_switchrhos1t.csv");
rhos2t = load("./results/results_kuramoto_switchrhos2t.csv");
rhos3t = load("./results/results_kuramoto_switchrhos3t.csv");
rhos4t = load("./results/results_kuramoto_switchrhos4t.csv");

spectrum_fft = load("./results/results_kuramoto_switchfft_spectrum.csv");
spectrum_fft_t = load("./results/results_kuramoto_switchfft_spectrum_t.csv");

%% Plotting
k = 10;

figure
subplot(211)
plot(1:250,spectrum1(k,:))
title("ISS")
grid on
subplot(212)
plot(1:250,spectrum1t(k,:))
grid on

figure
subplot(211)
plot(1:250,spectrum_fft(k,:))
title("FFT")
grid on
subplot(212)
plot(1:250,spectrum_fft_t(k,:))
grid on

%% check the rhos

figure
subplot(411)
plot(rhos1), hold on
plot(mean(rhos1).*ones(1,length(rhos1)),'r--')
grid on
subplot(412)
plot(rhos2), hold on
plot(mean(rhos2).*ones(1,length(rhos2)),'r--')
grid on
subplot(413)
plot(rhos3), hold on
plot(mean(rhos3).*ones(1,length(rhos3)),'r--')
grid on
subplot(414)
plot(rhos4), hold on
plot(mean(rhos4).*ones(1,length(rhos4)),'r--')
grid on

figure
subplot(411)
plot(rhos1t), hold on
plot(mean(rhos1t).*ones(1,length(rhos1t)),'r--')
grid on
subplot(412)
plot(rhos2t), hold on
plot(mean(rhos2t).*ones(1,length(rhos2t)),'r--')
grid on
subplot(413)
plot(rhos3t), hold on
plot(mean(rhos3t).*ones(1,length(rhos3t)),'r--')
grid on
subplot(414)
plot(rhos4t), hold on
plot(mean(rhos4t).*ones(1,length(rhos4t)),'r--')
grid on

%% Plot spectra

t1 = 1:size(spectrum1,1);
t2 = 1:size(spectrum1,2);

fs = 22;
lw2 = 1;

len = 90;

cax_max = 0.013;

% 
% figure('Units','normalized','Position',[.01 .01 .99 .99])
% 
% subplot(4,1,[1,2])
% surf(t1,t2,spectrum1')
% view(2)
% shading flat
% caxis([0 0.01])
% xlim([t1(1) t1(end)])
% ylim([2 t2(end)])
% title("Evolutionary inter spike spectrum of Kuramoto with switch")
% ylabel("Period [log_{10} a.u.]")
% set(gca,'Yscale','log')
% set(gca,'FontSize',fs)
% set(gca,'LineWidth',2)
% legend("\rho-thres = 0.8")
% grid on
% 
% subplot(4,1,[3,4])
% surf(t1,t2,spectrum2')
% view(2)
% shading flat
% caxis([0 0.01])
% xlim([t1(1) t1(end)])
% ylim([2 t2(end)])
% ylabel("Period [log_{10} a.u.]")
% xlabel("time [a.u.]")
% set(gca,'Yscale','log')
% set(gca,'FontSize',fs)
% set(gca,'LineWidth',2)
% legend("\rho-thres = 0.85")
% grid on
% 
% figure('Units','normalized','Position',[.01 .01 .99 .99])
% 
% subplot(4,1,[1,2])
% surf(t1,t2,spectrum3')
% view(2)
% shading flat
% caxis([0 0.01])
% xlim([t1(1) t1(end)])
% ylim([2 t2(end)])
% title("Evolutionary inter spike spectrum of Kuramoto with switch")
% ylabel("Period [log_{10} a.u.]")
% set(gca,'Yscale','log')
% set(gca,'FontSize',fs)
% set(gca,'LineWidth',2)
% legend("\rho-thres = 0.9")
% grid on
% 
% subplot(4,1,[3,4])
% surf(t1,t2,spectrum4')
% view(2)
% shading flat
% caxis([0 0.01])
% xlim([t1(1) t1(end)])
% ylim([2 t2(end)])
% ylabel("Period [log_{10} a.u.]")
% xlabel("time [a.u.]")
% set(gca,'Yscale','log')
% set(gca,'FontSize',fs)
% set(gca,'LineWidth',2)
% legend("\rho-thres = 0.95")
% grid on



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
legend("\rho-thres = 0.8")
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
legend("\rho-thres = 0.85")
grid on

figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(t1,t2,spectrum3')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
title("Evolutionary inter spike spectrogram of Kuramoto with switch")
ylabel("Period [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\rho-thres = 0.9")
grid on

subplot(4,1,[3,4])
surf(t1,t2,spectrum4')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
ylabel("Period [a.u.]")
xlabel("time [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\rho-thres = 0.95")
grid on


%

figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(t1,t2,spectrum1t')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
title("Evolutionary inter spike spectrogram of Kuramoto with switch TRANS")
ylabel("Period [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\rho-thres = 0.8")
grid on

subplot(4,1,[3,4])
surf(t1,t2,spectrum2t')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
ylabel("Period [a.u.]")
xlabel("time [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\rho-thres = 0.85")
grid on

figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(t1,t2,spectrum3t')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
title("Evolutionary inter spike spectrogram of Kuramoto with switch TRANS")
ylabel("Period [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\rho-thres = 0.9")
grid on

subplot(4,1,[3,4])
surf(t1,t2,spectrum4t')
view(2)
shading flat
caxis([0 0.01])
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
ylabel("Period [a.u.]")
xlabel("time [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
legend("\rho-thres = 0.95")
grid on

%

figure('Units','normalized','Position',[.01 .01 .99 .99])

subplot(4,1,[1,2])
surf(t1,t2,spectrum_fft')
view(2)
shading flat
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
title("Evolutionary FFT spectrogram of Kuramoto with switch")
ylabel("Period [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
grid on

subplot(4,1,[3,4])
surf(t1,t2,spectrum_fft_t')
view(2)
shading flat
xlim([t1(1) t1(end)])
ylim([2 t2(end)])
title("Evolutionary FFT spectrogram of Kuramoto with switch TRANS")
ylabel("Period [a.u.]")
xlabel("time [a.u.]")
set(gca,'FontSize',fs)
set(gca,'LineWidth',2)
grid on

%% Plot Wasserstein distance of all the spectra
% This has been computed in the script `compute_emd_system_with_transition.jl`

emd1 = load("./results/emd1.csv");
emd2 = load("./results/emd2.csv");
emd3 = load("./results/emd3.csv");
emd4 = load("./results/emd4.csv");

%%

imagesc(emd1<30), colormap([1 1 1; 0 0 0]), axis xy square
