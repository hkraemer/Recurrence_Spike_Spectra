%% Visualize results for Kuramoto with transition

% Computations have been carried out on the cluster (windowed spectrum determination)
% in the script `analyze_system_with_transition_cl.jl`

clear, clc
% load the results for the different rho-thresholds:

spectrum1 = load("./results/spectrum_0_8.csv");
spectrum2 = load("./results/spectrum_0_85.csv");
spectrum3 = load("./results/spectrum_0_9.csv");
spectrum4 = load("./results/spectrum_0_95.csv");

rhos1 = load("./results/rhos_0_8.csv");
rhos2 = load("./results/rhos_0_85.csv");
rhos3 = load("./results/rhos_0_9.csv");
rhos4 = load("./results/rhos_0_95.csv");

%% Plotting

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

%% Plot Wasserstein distance of all the spectra
% This has been computed in the script `compute_emd_system_with_transition.jl`

%% Load the EMD-Distance matrices
emd1 = load("./results/emd1.csv");
emd2 = load("./results/emd2.csv");
emd3 = load("./results/emd3.csv");
emd4 = load("./results/emd4.csv");

% emd1t = load("./results/emd1t.csv");
% emd2t = load("./results/emd2t.csv");
% emd3t = load("./results/emd3t.csv");
% emd4t = load("./results/emd4t.csv");
% 
% emd_fft = load("./results/emd_fft.csv");
% emd_fft_t = load("./results/emd_fft_t.csv");

%% save as .mat-files
save("./results/emd1.mat","emd1")
save("./results/emd2.mat","emd2")
save("./results/emd3.mat","emd3")
save("./results/emd4.mat","emd4")

% save("./results/emd1t.mat","emd1t")
% save("./results/emd2t.mat","emd2t")
% save("./results/emd3t.mat","emd3t")
% save("./results/emd4t.mat","emd4t")
% 
% save("./results/emd_fft.mat","emd_fft")
% save("./results/emd_fft_t.mat","emd_fft_t")

%% Load .mat files
emd1 = load("./results/emd1.mat");
emd2 = load("./results/emd2.mat");
emd3 = load("./results/emd3.mat");
emd4 = load("./results/emd4.mat");

emd1 = emd1.emd1;
emd2 = emd2.emd2;
emd3 = emd3.emd3;
emd4 = emd4.emd4;

% emd1t = load("./results/emd1t.mat");
% emd2t = load("./results/emd2t.mat");
% emd3t = load("./results/emd3t.mat");
% emd4t = load("./results/emd4t.mat");

% emd1t = emd1t.emd1t;
% emd2t = emd2t.emd2t;
% emd3t = emd3t.emd3t;
% emd4t = emd4t.emd4t;

% emd_fft = load("./results/emd_fft.mat");
% emd_fft_t = load("./results/emd_fft_t.mat");

% emd_fft = emd_fft.emd_fft;
% emd_fft_t = emd_fft_t.emd_fft_t;

%% Compute and plot RPs of EMDs

method = "var";
e  = 0.25;
RP = compute_rp_from_distance_matrix(emd1, method, e);

imagesc(RP), colormap([1 1 1; 0 0 0]), axis xy square




%%

%% Helper function

function y = compute_rp_from_distance_matrix(P, method, e)

len = length(P);
% compute recurrence matrix
if strcmp(method,"var")
    epsilon = quantile(P(:),e);
    y = double(P < epsilon); 
else
    q = quantile(P,e); % distance that corresponds to the fraction e of rec. points per column
    thresholds = repmat(q,len,1); % q has to be applied for each row in d
    % apply individual thresholds
    % apply threshold(s)
    y=double(P<thresholds);
end
end
